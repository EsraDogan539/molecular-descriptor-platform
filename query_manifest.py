"""Reproducible database-query manifest helpers for ChalMolDB."""

import hashlib
import json
from datetime import datetime, timezone

from release_metadata import DATABASE_VERSION, SCIENTIFIC_CORE_VERSION


QUERY_MANIFEST_SCHEMA_VERSION = "1.0"


def _canonical_json(value):
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )


def sha256_json(value):
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def result_set_sha256(record_ids):
    normalized = [str(value) for value in record_ids]
    return sha256_json(normalized)


def build_query_manifest(
    *,
    text_query,
    structure_query,
    standardized_structure_query,
    structure_mode,
    minimum_similarity,
    maximum_results,
    fingerprint_method,
    fingerprint_radius,
    fingerprint_bits,
    similarity_metric,
    filters,
    result_record_ids,
):
    """Build a machine-readable manifest for one database query state."""
    structure = {
        "input_smiles": structure_query or None,
        "standardized_smiles": standardized_structure_query or None,
        "mode": structure_mode if structure_query else None,
    }

    fingerprint = None
    if structure_query and structure_mode == "Similarity":
        fingerprint = {
            "method": fingerprint_method,
            "radius": int(fingerprint_radius),
            "bits": int(fingerprint_bits),
            "metric": similarity_metric,
            "minimum_similarity": float(minimum_similarity),
            "maximum_results": int(maximum_results),
        }

    query_configuration = {
        "text_query": text_query or None,
        "structure": structure,
        "fingerprint": fingerprint,
        "filters": filters,
        "database_version": DATABASE_VERSION,
        "scientific_core_version": SCIENTIFIC_CORE_VERSION,
    }

    record_ids = [str(value) for value in result_record_ids]

    return {
        "schema_version": QUERY_MANIFEST_SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "database_version": DATABASE_VERSION,
        "scientific_core_version": SCIENTIFIC_CORE_VERSION,
        "query": query_configuration,
        "result_count": len(record_ids),
        "query_checksum_sha256": sha256_json(query_configuration),
        "result_record_ids_sha256": result_set_sha256(record_ids),
    }


def query_manifest_json(manifest):
    return json.dumps(
        manifest,
        ensure_ascii=False,
        indent=2,
        sort_keys=True,
        default=str,
    )


def load_query_manifest(payload):
    """Parse and validate a ChalMolDB query-manifest payload.

    Accepts UTF-8 bytes, JSON text, or an already-decoded mapping.
    Validation is intentionally limited to fields needed for replay so that
    compatible manifests from earlier scientific-core releases remain usable.
    """
    if isinstance(payload, bytes):
        payload = payload.decode("utf-8")
    if isinstance(payload, str):
        try:
            manifest = json.loads(payload)
        except json.JSONDecodeError as exc:
            raise ValueError("The uploaded file is not valid JSON.") from exc
    elif isinstance(payload, dict):
        manifest = payload
    else:
        raise ValueError("Unsupported query-manifest payload.")

    if not isinstance(manifest, dict):
        raise ValueError("The query manifest must contain a JSON object.")

    schema_version = str(manifest.get("schema_version", "")).strip()
    if schema_version != QUERY_MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported query-manifest schema: {schema_version or 'missing'}. "
            f"Expected {QUERY_MANIFEST_SCHEMA_VERSION}."
        )

    query = manifest.get("query")
    if not isinstance(query, dict):
        raise ValueError("The query manifest does not contain a valid query configuration.")

    structure = query.get("structure")
    if structure is not None and not isinstance(structure, dict):
        raise ValueError("The query manifest contains an invalid structure configuration.")

    filters = query.get("filters")
    if filters is not None and not isinstance(filters, dict):
        raise ValueError("The query manifest contains an invalid filters configuration.")

    return manifest


def replay_configuration(manifest):
    """Return the replayable query configuration from a validated manifest."""
    manifest = load_query_manifest(manifest)
    query = manifest["query"]
    structure = query.get("structure") or {}
    fingerprint = query.get("fingerprint") or {}
    filters = query.get("filters") or {}

    return {
        "text_query": query.get("text_query") or "",
        "structure_query": structure.get("input_smiles") or "",
        "structure_mode": structure.get("mode") or "Exact",
        "minimum_similarity": fingerprint.get("minimum_similarity", 0.40),
        "maximum_results": fingerprint.get("maximum_results", 50),
        "filters": filters,
        "manifest_database_version": manifest.get("database_version"),
        "manifest_scientific_core_version": manifest.get("scientific_core_version"),
    }
