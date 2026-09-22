import pandas as pd

from public_labels import (
    public_collection_label,
    public_record_id,
    sanitize_public_dataframe,
)


def test_public_record_ids_remove_person_linked_prefixes():
    assert public_record_id("EROL_0001") == "DEV_0001"
    assert public_record_id("HAKAN_0001") == "EXT_0001"


def test_public_collection_labels_are_neutral():
    assert public_collection_label("Development/Training", "Erol dataset") == "Development collection"
    assert public_collection_label("External Validation", "Hakan dataset") == "External literature collection"


def test_public_dataframe_does_not_expose_person_linked_labels():
    df = pd.DataFrame(
        {
            "Record_ID": ["EROL_0001", "HAKAN_0001"],
            "Dataset_Owner": ["Erol dataset", "Hakan dataset"],
            "Split_Role": ["Development/Training", "External Validation"],
        }
    )
    public = sanitize_public_dataframe(df)
    text = public.to_csv(index=False).lower()
    assert "erol" not in text
    assert "hakan" not in text
    assert "dev_0001" in text
    assert "ext_0001" in text
