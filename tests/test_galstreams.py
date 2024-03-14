def test_import():
    import galstreams as gs

    mws = gs.MWStreams()
    assert len(mws.all_track_names()) > 0
