"""Record label (TEXT) selection in binary layer files.

The files are written here with struct, independently of the flopy writers.
Each value encodes its label, layer, stress period, time step, node and
occurrence, so a record read from the wrong position or label is detected.
"""

import struct

import numpy as np
import pytest

from flopy.utils import HeadFile, HeadUFile, UcnFile

# HeadUFile layers (ilay, first node, last node), of different sizes
LAYERS = ((1, 1, 3), (2, 4, 8))
# (kper, kstp, pertim, totim); kstp and pertim restart in the second period
CALENDAR = ((1, 1, 5.0, 5.0), (1, 2, 10.0, 10.0), (2, 1, 5.0, 15.0), (2, 2, 10.0, 20.0))
LABELS = {1: b"          CONC01", 2: b"          CONC02", 0: b"           HEADU"}
FMT = {"single": "f", "double": "d"}
DTYPE = {"single": np.float32, "double": np.float64}


def value(label, ilay, kper, kstp, node, occurrence=0):
    return 10000 * occurrence + 1000 * label + 100 * ilay + 10 * kper + kstp + node / 8


def write_headu(path, blocks, precision):
    """Write an unstructured (ULASAVU) file and return the records written.

    blocks is a list of (label, occurrence) written for each layer and time
    step of CALENDAR.
    """
    f = FMT[precision]
    records = []
    with open(path, "wb") as fh:
        for kper, kstp, pertim, totim in CALENDAR:
            for label, occurrence in blocks:
                for ilay, nstrt, nend in LAYERS:
                    data = [
                        value(label, ilay, kper, kstp, node, occurrence)
                        for node in range(nstrt, nend + 1)
                    ]
                    fh.write(
                        struct.pack(
                            f"<2i2{f}16s3i",
                            kstp,
                            kper,
                            pertim,
                            totim,
                            LABELS[label],
                            nstrt,
                            nend,
                            ilay,
                        )
                    )
                    offset = fh.tell()
                    fh.write(struct.pack(f"<{len(data)}{f}", *data))
                    records.append(
                        {
                            "kper": kper,
                            "kstp": kstp,
                            "pertim": pertim,
                            "totim": totim,
                            "text": LABELS[label],
                            "label": label,
                            "occurrence": occurrence,
                            "ilay": ilay,
                            "nstrt": nstrt,
                            "nend": nend,
                            "offset": offset,
                            "data": np.array(data, dtype=DTYPE[precision]),
                        }
                    )
    return records


def check_index(hds, expected):
    assert len(hds) == len(expected)
    assert len(hds.iposarray) == len(hds.recordarray) == len(hds.headers)
    assert hds.recordarray.tolist() == [
        (
            r["kstp"],
            r["kper"],
            r["pertim"],
            r["totim"],
            r["text"],
            r["nstrt"],
            r["nend"],
            r["ilay"],
        )
        for r in expected
    ]
    assert hds.iposarray.tolist() == [r["offset"] for r in expected]
    assert hds.headers.index.tolist() == [r["offset"] for r in expected]
    assert hds.headers["text"].tolist() == [
        r["text"].decode().strip() for r in expected
    ]
    assert hds.nlay == 2
    assert hds.get_times() == [5.0, 10.0, 15.0, 20.0]
    assert hds.kstpkper == [(1, 1), (2, 1), (1, 2), (2, 2)]
    assert hds.get_kstpkper() == [(0, 0), (1, 0), (0, 1), (1, 1)]


def check_values(hds, expected, precision):
    """Each layer and time step returns the last expected record."""
    for (kper, kstp, pertim, totim), kstpkper in zip(CALENDAR, hds.get_kstpkper()):
        last = {}
        for r in expected:
            if r["totim"] == totim:
                last[r["ilay"]] = r["data"]
        for data in (hds.get_data(kstpkper=kstpkper), hds.get_data(totim=totim)):
            assert len(data) == 2
            for ilay in (1, 2):
                assert data[ilay - 1].dtype == DTYPE[precision]
                np.testing.assert_array_equal(data[ilay - 1], last[ilay])


@pytest.mark.parametrize("precision", ["single", "double"])
@pytest.mark.parametrize(
    "label, text",
    [(0, "headu"), (1, "conc01"), (1, "conc")],
    ids=["headu", "conc01", "conc"],
)
def test_headufile_single_label(tmp_path, precision, label, text):
    records = write_headu(tmp_path / "single.bin", [(label, 0)], precision)
    hds = HeadUFile(tmp_path / "single.bin", text=text, precision=precision)
    check_index(hds, records)
    check_values(hds, records, precision)
    # ragged layers: the header holds the first and last node of each layer
    assert [len(d) for d in hds.get_data()] == [3, 5]
    hds.close()


@pytest.mark.parametrize("precision", ["single", "double"])
@pytest.mark.parametrize("label", [1, 2])
def test_headufile_interleaved_labels(tmp_path, precision, label):
    records = write_headu(tmp_path / "conc.bin", [(1, 0), (2, 0)], precision)
    expected = [r for r in records if r["label"] == label]
    hds = HeadUFile(tmp_path / "conc.bin", text=f"conc0{label}", precision=precision)
    check_index(hds, expected)
    check_values(hds, expected, precision)
    # metadata come from the selected records, also when the first record
    # in the file has another label
    assert (hds.ncol, hds.nrow) == (1, 3)
    other = 3 - label
    for data in (hds.get_data(totim=t) for t in hds.get_times()):
        for d in data:
            assert np.all(d // 1000 % 10 == label)
            assert not np.any(d // 1000 % 10 == other)
    # time series and all data of one layer for the selected label
    node = 5  # zero-based, second node of layer 2
    ts = hds.get_ts(idx=node)
    np.testing.assert_array_equal(ts[:, 0], [5.0, 10.0, 15.0, 20.0])
    np.testing.assert_array_equal(
        ts[:, 1],
        [value(label, 2, kper, kstp, node + 1) for kper, kstp, _, _ in CALENDAR],
    )
    alldata = hds.get_alldata(mflay=1)
    np.testing.assert_array_equal(
        alldata, [r["data"] for r in expected if r["ilay"] == 2]
    )
    hds.close()


@pytest.mark.parametrize("precision", ["single", "double"])
def test_headufile_generic_text_with_several_labels_raises(tmp_path, precision):
    write_headu(tmp_path / "conc.bin", [(1, 0), (2, 0)], precision)
    hds = HeadUFile(tmp_path / "conc.bin", text="conc", precision=precision)
    # all records match the generic text and are indexed
    assert len(hds) == 16
    assert hds.get_times() == [5.0, 10.0, 15.0, 20.0]
    match = "more specific text"
    with pytest.raises(ValueError, match=match):
        hds.get_data(kstpkper=(0, 0))
    with pytest.raises(ValueError, match=match):
        hds.get_data(totim=20.0)
    with pytest.raises(ValueError, match=match):
        hds.get_data()
    with pytest.raises(ValueError, match=match):
        hds.get_alldata(mflay=0)
    with pytest.raises(ValueError, match=match):
        hds.get_ts(idx=0)
    hds.close()


@pytest.mark.parametrize("precision", ["single", "double"])
def test_headufile_repeated_label_keeps_last_occurrence(tmp_path, precision):
    # CONC01 is written twice per time step, with CONC02 in between
    records = write_headu(tmp_path / "conc.bin", [(1, 0), (2, 0), (1, 1)], precision)

    hds = HeadUFile(tmp_path / "conc.bin", text="conc01", precision=precision)
    expected = [r for r in records if r["label"] == 1]
    check_index(hds, expected)
    # the last occurrence of the label is returned for each layer
    check_values(hds, expected, precision)
    for data in (hds.get_data(totim=t) for t in hds.get_times()):
        for d in data:
            assert np.all(d // 10000 == 1)
            assert np.all(d // 1000 % 10 == 1)
    hds.close()

    hds = HeadUFile(tmp_path / "conc.bin", text="conc02", precision=precision)
    expected = [r for r in records if r["label"] == 2]
    check_index(hds, expected)
    check_values(hds, expected, precision)
    hds.close()


def test_headufile_text_not_found_raises(tmp_path):
    write_headu(tmp_path / "conc.bin", [(1, 0), (2, 0)], "single")
    with pytest.raises(ValueError, match="not found"):
        HeadUFile(tmp_path / "conc.bin", text="headu", precision="single")


@pytest.mark.parametrize("precision", ["single", "double"])
def test_headufile_truncated_file_raises(tmp_path, precision):
    records = write_headu(tmp_path / "conc.bin", [(1, 0), (2, 0)], precision)
    data = (tmp_path / "conc.bin").read_bytes()

    # last data array cut short: the index is built, the last array fails
    (tmp_path / "short_data.bin").write_bytes(data[:-3])
    hds = HeadUFile(tmp_path / "short_data.bin", text="conc02", precision=precision)
    assert hds.get_times() == [5.0, 10.0, 15.0, 20.0]
    (layer2,) = [
        r["data"]
        for r in records
        if (r["label"], r["ilay"], r["totim"]) == (2, 2, 15.0)
    ]
    np.testing.assert_array_equal(hds.get_data(totim=15.0)[1], layer2)
    with pytest.raises(EOFError):
        hds.get_data(totim=20.0)
    hds.close()

    # last header cut short
    (tmp_path / "short_header.bin").write_bytes(data[: records[-1]["offset"] - 10])
    with pytest.raises(EOFError):
        HeadUFile(tmp_path / "short_header.bin", text="conc01", precision=precision)


def write_structured(path, labels, precision, ucn=False):
    """Write a structured head (or MT3D concentration) file, 2 layers of 2 x 3."""
    f = FMT[precision]
    nrow, ncol = 2, 3
    records = []
    with open(path, "wb") as fh:
        for kper, kstp, pertim, totim in CALENDAR:
            for label, text in labels:
                for ilay in (1, 2):
                    data = [
                        value(label, ilay, kper, kstp, node)
                        for node in range(1, nrow * ncol + 1)
                    ]
                    if ucn:
                        hdr = struct.pack(
                            f"<3i{f}16s3i", 7, kstp, kper, totim, text, ncol, nrow, ilay
                        )
                    else:
                        hdr = struct.pack(
                            f"<2i2{f}16s3i",
                            kstp,
                            kper,
                            pertim,
                            totim,
                            text,
                            ncol,
                            nrow,
                            ilay,
                        )
                    fh.write(hdr)
                    offset = fh.tell()
                    fh.write(struct.pack(f"<{len(data)}{f}", *data))
                    records.append(
                        {
                            "totim": totim,
                            "text": text,
                            "ilay": ilay,
                            "offset": offset,
                            "data": np.array(data, dtype=DTYPE[precision]).reshape(
                                nrow, ncol
                            ),
                        }
                    )
    return records


def check_structured(obj, expected, precision):
    assert obj.iposarray.tolist() == [r["offset"] for r in expected]
    assert obj.recordarray["text"].tolist() == [r["text"] for r in expected]
    assert len(obj.headers) == len(expected)
    assert (obj.nlay, obj.nrow, obj.ncol) == (2, 2, 3)
    assert obj.get_times() == [5.0, 10.0, 15.0, 20.0]
    for totim in obj.get_times():
        data = obj.get_data(totim=totim)
        assert data.shape == (2, 2, 3)
        assert data.dtype == DTYPE[precision]
        for r in expected:
            if r["totim"] == totim:
                np.testing.assert_array_equal(data[r["ilay"] - 1], r["data"])


@pytest.mark.parametrize("precision", ["single", "double"])
def test_headfile_single_label(tmp_path, precision):
    records = write_structured(
        tmp_path / "model.hds", [(0, b"            HEAD")], precision
    )
    with HeadFile(tmp_path / "model.hds", precision=precision) as hds:
        check_structured(hds, records, precision)


@pytest.mark.parametrize("precision", ["single", "double"])
@pytest.mark.parametrize("text", ["head", "drawdown"])
def test_headfile_interleaved_labels(tmp_path, precision, text):
    labels = [(1, b"            HEAD"), (2, b"        DRAWDOWN")]
    records = write_structured(tmp_path / "model.bin", labels, precision)
    expected = [r for r in records if text.upper().encode() in r["text"]]
    with HeadFile(tmp_path / "model.bin", text=text, precision=precision) as hds:
        check_structured(hds, expected, precision)


@pytest.mark.parametrize("precision", ["single", "double"])
def test_ucnfile_single_label(tmp_path, precision):
    records = write_structured(
        tmp_path / "MT3D001.UCN", [(1, b"CONCENTRATION   ")], precision, ucn=True
    )
    with UcnFile(tmp_path / "MT3D001.UCN", precision=precision) as ucn:
        check_structured(ucn, records, precision)
        assert ucn.recordarray["ntrans"].tolist() == [7] * len(records)
