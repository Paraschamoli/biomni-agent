# biomni_agent/biomedical_tools.py
"""Lightweight biomedical tools for the Biomni agent (Production Safe Version)."""

import re
import tempfile
import asyncio
import base64
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Union

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import scanpy as sc
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

from agno.tools import tool


# REQUIRED by NCBI - type ignore for Entrez.email (BioPython typing issue)
Entrez.email = "paraschamoli2592002@gmail.com"  # type: ignore[assignment]


# ==============================
# SEQUENCE VALIDATION
# ==============================

def validate_sequence(sequence: str, seq_type: str = "dna") -> Dict[str, Any]:
    sequence = sequence.upper().strip()

    if not sequence:
        return {"valid": False, "error": "Empty sequence"}

    if len(sequence) < 20:
        return {"valid": False, "error": "Sequence must be at least 20 characters"}

    if seq_type == "dna":
        valid_chars = set("ACGTN")
    elif seq_type == "rna":
        valid_chars = set("ACGUN")
    elif seq_type == "protein":
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY*")
    else:
        return {"valid": False, "error": "Invalid sequence type"}

    invalid = set(sequence) - valid_chars
    if invalid:
        return {"valid": False, "error": f"Invalid characters: {''.join(invalid)}"}

    return {"valid": True, "sequence": sequence}


def validate_smiles(smiles: str) -> Dict[str, Any]:
    smiles = smiles.strip()

    if not smiles:
        return {"valid": False, "error": "Empty SMILES string"}

    pattern = r'^[A-Za-z0-9@+\-\[\]\(\)\\\/%=#$.~]*$'
    if not re.match(pattern, smiles):
        return {"valid": False, "error": "Invalid SMILES format"}

    return {"valid": True, "smiles": smiles}


# ==============================
# BLAST TOOL
# ==============================

@tool
async def run_blast_sequence_async(
    sequence: str,
    program: str = "blastn",
    database: str = "nr"
) -> Dict[str, Any]:

    try:
        seq_type = "dna" if program != "blastp" else "protein"
        validation = validate_sequence(sequence, seq_type)
        if not validation["valid"]:
            return {"success": False, "error": validation["error"]}

        async with asyncio.timeout(180):  # 3 min timeout
            result_handle = await asyncio.to_thread(
                NCBIWWW.qblast,
                program,
                database,
                validation["sequence"],
                format_type="XML"
            )

        blast_records = NCBIXML.parse(result_handle)

        hits = []
        for record in blast_records:
            for alignment in record.alignments[:5]:
                for hsp in alignment.hsps[:1]:
                    hits.append({
                        "title": alignment.hit_def,
                        "e_value": float(hsp.expect),
                        "identity_percent": round(
                            (hsp.identities / hsp.align_length) * 100, 2
                        ),
                        "score": float(hsp.score)
                    })

        result_handle.close()

        return {
            "success": True,
            "query_length": len(validation["sequence"]),
            "num_hits": len(hits),
            "hits": hits,
            "summary": f"Found {len(hits)} BLAST hits."
        }

    except asyncio.TimeoutError:
        return {"success": False, "error": "BLAST search timed out."}
    except Exception as e:
        return {"success": False, "error": f"BLAST error: {str(e)}"}


# ==============================
# scRNA-seq TOOL
# ==============================

@tool
async def analyze_scrnaseq_async(
    data_path: Union[str, Path],
    min_genes: int = 200,
    min_cells: int = 3,
    n_pcs: int = 50,
    n_neighbors: int = 10
) -> Dict[str, Any]:

    try:
        # Convert to Path object for file operations
        path_obj = Path(data_path)
        if not path_obj.exists():
            return {"success": False, "error": "File not found"}

        adata = sc.read(str(path_obj))  # Convert to string for scanpy

        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

        sc.pp.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1))
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)
        sc.tl.leiden(adata)

        # Safely get cluster counts - type ignore for scanpy typing issue
        if "leiden" in adata.obs:
            clusters = adata.obs["leiden"].value_counts().to_dict()  # type: ignore[union-attr]
        else:
            clusters = {}

        return {
            "success": True,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_clusters": len(clusters),
            "cluster_sizes": clusters,
            "summary": f"Identified {len(clusters)} clusters."
        }

    except Exception as e:
        return {"success": False, "error": f"scRNA-seq error: {str(e)}"}


# ==============================
# ADMET TOOL (LAZY LOAD FIX)
# ==============================

@tool
async def predict_admet_properties_async(
    smiles: str,
    molecule_name: Optional[str] = None
) -> Dict[str, Any]:

    try:
        validation = validate_smiles(smiles)
        if not validation["valid"]:
            return {"success": False, "error": validation["error"]}

        # 🔥 Lazy import (prevents startup crash)
        try:
            from chembl_webresource_client.new_client import new_client
        except Exception as e:
            return {
                "success": False,
                "error": f"ChEMBL client unavailable: {str(e)}"
            }

        molecule = new_client.molecule
        activity = new_client.activity

        async with asyncio.timeout(60):
            molecules = await asyncio.to_thread(
                lambda: list(
                    molecule.filter(
                        molecule_structures__canonical_smiles__flexmatch=validation["smiles"]
                    )
                )
            )

        if not molecules:
            return {"success": True, "message": "Molecule not found"}

        mol = molecules[0]

        async with asyncio.timeout(30):
            activities = await asyncio.to_thread(
                lambda: list(
                    activity.filter(
                        molecule_chembl_id=mol["molecule_chembl_id"],
                        limit=20
                    )
                )
            )

        return {
            "success": True,
            "molecule": mol.get("pref_name", "Unknown"),
            "chembl_id": mol["molecule_chembl_id"],
            "activities_found": len(activities)
        }

    except asyncio.TimeoutError:
        return {"success": False, "error": "ChEMBL query timed out"}
    except Exception as e:
        return {"success": False, "error": f"ADMET error: {str(e)}"}


# ==============================
# PDF TOOL
# ==============================

@tool
async def generate_pdf_report_async(
    title: str,
    content: str,
    author: str = "BioOmni Agent"
) -> Dict[str, Any]:

    try:
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
            pdf_path = tmp.name

        c = canvas.Canvas(pdf_path, pagesize=letter)
        width, height = letter

        c.setTitle(title)
        c.setAuthor(author)

        y = height - 50
        c.setFont("Helvetica-Bold", 16)
        c.drawString(50, y, title)

        y -= 30
        c.setFont("Helvetica", 10)
        c.drawString(50, y, f"Generated: {datetime.now()}")

        y -= 40
        c.setFont("Helvetica", 12)

        for line in content.split("\n"):
            if y < 50:
                c.showPage()
                c.setFont("Helvetica", 12)
                y = height - 50
            c.drawString(50, y, line[:100])
            y -= 18

        c.save()

        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()

        Path(pdf_path).unlink(missing_ok=True)

        return {
            "success": True,
            "size_bytes": len(pdf_bytes),
            "summary": "PDF generated successfully"
        }

    except Exception as e:
        return {"success": False, "error": f"PDF error: {str(e)}"}