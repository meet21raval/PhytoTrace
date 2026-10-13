
import csv
import io
import os
import re
import sys
import time
import queue
import threading
import traceback
from pathlib import Path
from urllib.parse import quote, urljoin

import requests
from bs4 import BeautifulSoup
import pandas as pd

# RDKit is used only for descriptors computed from a retrieved SMILES.
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, QED
    RDKIT_OK = True
except Exception:
    RDKIT_OK = False

import tkinter as tk
from tkinter import ttk, filedialog, messagebox


APP_TITLE = "Plant Phytoconstituent Retriever — Version 1"
IMPPAT_BASE = "https://cb.imsc.res.in/imppat"
IMPPAT_PLANT_URL = IMPPAT_BASE + "/phytochemical/{plant}"
IMPPAT_DETAIL_URL = IMPPAT_BASE + "/phytochemical-detailedpage/{phy_id}"

NPASS_DOWNLOAD_URL = "https://bidd.group/NPASS/downloadnpass.html"

HEADERS = {
    "User-Agent": "Plant-Phytoconstituent-Retriever/1.0 (research data retrieval)",
    "Accept-Language": "en-US,en;q=0.9",
}

SESSION = requests.Session()
SESSION.headers.update(HEADERS)


def clean_text(x):
    if x is None:
        return ""
    return re.sub(r"\s+", " ", str(x)).strip()


def norm_key(x):
    return re.sub(r"[^a-z0-9]+", "", clean_text(x).lower())


def safe_filename(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "plant"


def get(url, timeout=40):
    r = SESSION.get(url, timeout=timeout)
    r.raise_for_status()
    return r


def parse_external_ids(text):
    out = {}
    for token in re.split(r"[,;]", clean_text(text)):
        token = token.strip()
        if ":" in token:
            k, v = token.split(":", 1)
            out[clean_text(k)] = clean_text(v)
    return out


def parse_imppat_detail(phy_id):
    url = IMPPAT_DETAIL_URL.format(phy_id=quote(phy_id, safe=""))
    r = get(url)
    soup = BeautifulSoup(r.text, "html.parser")
    text = soup.get_text("\n", strip=True)

    data = {
        "IMPPAT_ID": phy_id,
        "IMPPAT_detail_URL": url,
        "Phytochemical_name": "",
        "Synonyms": "",
        "External_identifiers": "",
        "PubChem_CID": "",
        "SMILES": "",
        "InChI": "",
        "InChIKey": "",
        "DeepSMILES": "",
        "Functional_groups": "",
        "ClassyFire_Kingdom": "",
        "ClassyFire_Superclass": "",
        "ClassyFire_Class": "",
        "ClassyFire_Subclass": "",
        "NPClassifier_Pathway": "",
        "NPClassifier_Superclass": "",
        "NPClassifier_Class": "",
        "NP_Likeness": "",
    }

    # Current IMPPAT 3.0 pages have a predictable text layout.
    patterns = {
        "Phytochemical_name": r"Phytochemical name:\s*(.+)",
        "Synonyms": r"Synonymous chemical names:\s*(.+?)(?=\nExternal chemical identifiers:)",
        "External_identifiers": r"External chemical identifiers:\s*(.+?)(?=\nChemical structure information)",
        "SMILES": r"SMILES:\s*(.+)",
        "InChI": r"InChI:\s*(.+)",
        "InChIKey": r"InChIKey:\s*(.+)",
        "DeepSMILES": r"DeepSMILES:\s*(.+)",
        "Functional_groups": r"Functional groups:\s*(.+?)(?=\nMolecular scaffolds)",
        "ClassyFire_Kingdom": r"ClassyFire Kingdom:\s*(.+)",
        "ClassyFire_Superclass": r"ClassyFire Superclass:\s*(.+?)\s*ClassyFire Class:",
        "ClassyFire_Class": r"ClassyFire Class:\s*(.+)",
        "ClassyFire_Subclass": r"ClassyFire Subclass:\s*(.+)",
        "NPClassifier_Pathway": r"NP Classifier Biosynthetic pathway:\s*(.+)",
        "NPClassifier_Superclass": r"NP Classifier Superclass:\s*(.+)",
        "NPClassifier_Class": r"NP Classifier Class:\s*(.+)",
        "NP_Likeness": r"NP-Likeness score:\s*([0-9.+-]+)",
    }

    for key, pat in patterns.items():
        m = re.search(pat, text, flags=re.I)
        if m:
            data[key] = clean_text(m.group(1))

    ext = parse_external_ids(data["External_identifiers"])
    for k, v in ext.items():
        if norm_key(k) in {"cid", "pubchemcid"}:
            data["PubChem_CID"] = v

    # More reliable extraction of SMILES/InChI blocks from the page text.
    for label, key in [("SMILES:", "SMILES"), ("InChI:", "InChI"),
                       ("InChIKey:", "InChIKey"), ("DeepSMILES:", "DeepSMILES")]:
        if not data[key]:
            pos = text.find(label)
            if pos >= 0:
                tail = text[pos + len(label):].splitlines()
                if tail:
                    data[key] = clean_text(tail[0])

    return data


def compute_rdkit(smiles):
    result = {
        "RDKit_valid_structure": "",
        "Molecular_formula": "",
        "MW_g_mol": "",
        "Exact_mass": "",
        "XlogP": "",
        "TPSA_A2": "",
        "HBD": "",
        "HBA": "",
        "Rotatable_bonds": "",
        "Heavy_atoms": "",
        "Carbon_atoms": "",
        "Nitrogen_atoms": "",
        "Oxygen_atoms": "",
        "Sulfur_atoms": "",
        "Halogen_atoms": "",
        "Ring_count": "",
        "Aromatic_ring_count": "",
        "Fraction_CSP3": "",
        "QED": "",
    }
    if not RDKIT_OK or not smiles:
        return result

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["RDKit_valid_structure"] = "No"
            return result

        result["RDKit_valid_structure"] = "Yes"
        result["Molecular_formula"] = rdMolDescriptors.CalcMolFormula(mol)
        result["MW_g_mol"] = round(Descriptors.MolWt(mol), 4)
        result["Exact_mass"] = round(Descriptors.ExactMolWt(mol), 6)
        result["XlogP"] = round(Crippen.MolLogP(mol), 4)
        result["TPSA_A2"] = round(rdMolDescriptors.CalcTPSA(mol), 4)
        result["HBD"] = Lipinski.NumHDonors(mol)
        result["HBA"] = Lipinski.NumHAcceptors(mol)
        result["Rotatable_bonds"] = Lipinski.NumRotatableBonds(mol)
        result["Heavy_atoms"] = Lipinski.HeavyAtomCount(mol)
        result["Carbon_atoms"] = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
        result["Nitrogen_atoms"] = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")
        result["Oxygen_atoms"] = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "O")
        result["Sulfur_atoms"] = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "S")
        result["Halogen_atoms"] = sum(
            1 for a in mol.GetAtoms() if a.GetSymbol() in {"F", "Cl", "Br", "I"}
        )
        result["Ring_count"] = Lipinski.RingCount(mol)
        result["Aromatic_ring_count"] = rdMolDescriptors.CalcNumAromaticRings(mol)
        result["Fraction_CSP3"] = round(rdMolDescriptors.CalcFractionCSP3(mol), 4)
        result["QED"] = round(QED.qed(mol), 4)
    except Exception:
        result["RDKit_valid_structure"] = "Error"
    return result


def extract_imppat_plant(plant_name, log=lambda x: None):
    plant_url = IMPPAT_PLANT_URL.format(plant=quote(plant_name.strip(), safe=""))
    log(f"IMPPAT: opening {plant_url}")
    r = get(plant_url)
    soup = BeautifulSoup(r.text, "html.parser")

    title = soup.find("h3")
    resolved_name = clean_text(title.get_text(" ", strip=True)) if title else plant_name

    # Plant metadata from the top section.
    page_text = soup.get_text("\n", strip=True)
    family = ""
    kingdom = ""
    group = ""
    synonyms = ""
    system_medicine = ""
    for label, target in [
        ("Kingdom:", "kingdom"), ("Family:", "family"), ("Group:", "group"),
        ("Synonymous names:", "synonyms"), ("System of medicine:", "system_medicine")
    ]:
        m = re.search(re.escape(label) + r"\s*(.+)", page_text)
        if m:
            locals()[target] = clean_text(m.group(1))

    records = []
    seen = set()

    # Table rows are the authoritative plant-part-phytochemical associations.
    for tr in soup.find_all("tr"):
        cells = tr.find_all(["td", "th"])
        if len(cells) < 4:
            continue
        row_text = [clean_text(c.get_text(" ", strip=True)) for c in cells]
        links = tr.find_all("a", href=True)

        phy_link = None
        phy_id = None
        phy_name = ""
        for a in links:
            txt = clean_text(a.get_text(" ", strip=True))
            href = a.get("href", "")
            if re.search(r"IMPP?HY\d+|IMPPAT3_PHYID\d+", txt, re.I) or \
               re.search(r"phytochemical-detailedpage", href, re.I):
                phy_link = a
                phy_id = txt if re.match(r"IMPPAT3_PHYID\d+|IMPHY\d+", txt, re.I) else None
                if not phy_id:
                    m = re.search(r"(IMPPAT3_PHYID\d+|IMPHY\d+)", href, re.I)
                    if m:
                        phy_id = m.group(1)
                break

        if not phy_id or not phy_link:
            continue

        # Determine plant part and name from table position.
        plant_part = row_text[1] if len(row_text) >= 2 else ""
        name_link = phy_link
        if not phy_name:
            phy_name = clean_text(name_link.get_text(" ", strip=True))

        key = (plant_part.lower(), phy_id.lower())
        if key in seen:
            continue
        seen.add(key)

        ref_text = row_text[4] if len(row_text) >= 5 else ""

        records.append({
            "Plant": resolved_name,
            "Plant_part": plant_part,
            "IMPPAT_ID": phy_id,
            "Phytochemical_name": phy_name,
            "IMPPAT_reference": ref_text,
            "IMPPAT_plant_URL": plant_url,
            "Plant_family": family,
            "Plant_kingdom": kingdom,
            "Plant_group": group,
            "Plant_synonyms": synonyms,
            "System_of_medicine": system_medicine,
        })

    log(f"IMPPAT: found {len(records)} plant-part/compound association rows")

    # Fetch unique compound detail pages in parallel.
    unique_ids = list(dict.fromkeys(x["IMPPAT_ID"] for x in records))
    detail_map = {}
    from concurrent.futures import ThreadPoolExecutor, as_completed
    with ThreadPoolExecutor(max_workers=6) as ex:
        futures = {ex.submit(parse_imppat_detail, pid): pid for pid in unique_ids}
        done = 0
        for fut in as_completed(futures):
            pid = futures[fut]
            done += 1
            try:
                detail_map[pid] = fut.result()
            except Exception as e:
                detail_map[pid] = {"IMPPAT_ID": pid, "IMPPAT_detail_URL": IMPPAT_DETAIL_URL.format(phy_id=pid)}
                log(f"IMPPAT detail failed for {pid}: {e}")
            if done % 10 == 0 or done == len(futures):
                log(f"IMPPAT: compound details {done}/{len(futures)}")

    output = []
    for rec in records:
        d = detail_map.get(rec["IMPPAT_ID"], {})
        merged = dict(rec)
        merged.update(d)
        merged.update(compute_rdkit(merged.get("SMILES", "")))
        merged["Source"] = "IMPPAT"
        output.append(merged)

    return pd.DataFrame(output)


def find_download_links(html):
    soup = BeautifulSoup(html, "html.parser")
    items = []
    for tr in soup.find_all("tr"):
        txt = clean_text(tr.get_text(" ", strip=True))
        links = [urljoin(NPASS_DOWNLOAD_URL, a["href"]) for a in tr.find_all("a", href=True)]
        if links:
            items.append((txt, links))
    return items


def download_npass_required_files(cache_dir, log=lambda x: None):
    """
    NPASS provides downloadable 2026 datasets. This function discovers the current
    links from the official download page rather than hard-coding changing URLs.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    r = get(NPASS_DOWNLOAD_URL)
    rows = find_download_links(r.text)

    wanted = {}
    for text, links in rows:
        low = text.lower()
        if "natural products-species source" in low or "species source" in low:
            wanted.setdefault("species_source", links[0])
        if "natural products-general information" in low or "natural products-general" in low:
            wanted.setdefault("general", links[0])
        if "natural products-structure" in low:
            wanted.setdefault("structure", links[0])

    if not wanted:
        raise RuntimeError(
            "Could not discover NPASS download links automatically. "
            "Use the NPASS website's 2026 download page and place the downloaded "
            "General Information, Structure, and Species Source files in the "
            "NPASS_cache folder."
        )

    paths = {}
    for key, url in wanted.items():
        ext = ".zip" if ".zip" in url.lower() else ".dat"
        path = cache_dir / f"npass_{key}{ext}"
        if not path.exists():
            log(f"NPASS: downloading {key} dataset...")
            rr = get(url, timeout=180)
            path.write_bytes(rr.content)
        paths[key] = path
    return paths


def read_table_file(path):
    """
    Flexible reader for NPASS text/CSV/TSV/ZIP files.
    Returns a DataFrame or None.
    """
    import zipfile
    candidates = []
    p = Path(path)

    if p.suffix.lower() == ".zip":
        with zipfile.ZipFile(p) as z:
            for n in z.namelist():
                if n.lower().endswith((".txt", ".tsv", ".csv", ".tab")):
                    candidates.append((n, z.read(n)))
    else:
        candidates.append((p.name, p.read_bytes()))

    for name, raw in candidates:
        for encoding in ("utf-8", "utf-8-sig", "latin-1"):
            try:
                txt = raw.decode(encoding)
                sample = txt[:5000]
                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters="\t,;|")
                    sep = dialect.delimiter
                except Exception:
                    sep = "\t"
                df = pd.read_csv(io.StringIO(txt), sep=sep, dtype=str, low_memory=False)
                if len(df.columns) >= 2:
                    return df
            except Exception:
                continue
    return None


def choose_column(df, patterns):
    for col in df.columns:
        nk = norm_key(col)
        for pat in patterns:
            if pat in nk:
                return col
    return None


def extract_npass_from_cache(plant_name, cache_dir, log=lambda x: None):
    paths = download_npass_required_files(cache_dir, log=log)

    species = read_table_file(paths["species_source"])
    general = read_table_file(paths["general"])
    structure = read_table_file(paths["structure"])

    if species is None:
        raise RuntimeError("NPASS species-source dataset could not be parsed.")

    name_col = choose_column(species, ["organismname", "speciesname", "organism", "species"])
    npid_col = choose_column(species, ["npid", "naturalproductid", "npcid", "compoundid"])

    if not name_col or not npid_col:
        raise RuntimeError(
            "NPASS species-source file format was not recognized. "
            "The cache files may need a one-time column mapping update."
        )

    target = plant_name.strip().lower()
    mask = species[name_col].fillna("").str.strip().str.lower().eq(target)
    hits = species.loc[mask].copy()

    if hits.empty:
        # Second pass: normalized whitespace/case matching.
        normalized_target = re.sub(r"\s+", " ", target)
        vals = species[name_col].fillna("").map(lambda x: re.sub(r"\s+", " ", str(x).strip().lower()))
        hits = species.loc[vals.eq(normalized_target)].copy()

    if hits.empty:
        return pd.DataFrame(columns=[
            "Plant", "NPASS_ID", "NPASS_name", "NPASS_SMILES",
            "NPASS_Molecular_Weight", "Source"
        ])

    ids = set(hits[npid_col].dropna().astype(str).str.strip())

    # Identify IDs/names/SMILES in general and structure tables.
    def prep(df):
        if df is None:
            return None, None
        idc = choose_column(df, ["npid", "naturalproductid", "npcid", "compoundid"])
        namec = choose_column(df, ["commonname", "compoundname", "naturalproductname", "name"])
        smic = choose_column(df, ["canonicalsmiles", "smiles"])
        mwc = choose_column(df, ["molecularweight", "mw"])
        return idc, (namec, smic, mwc)

    gid, gcols = prep(general)
    sid, scols = prep(structure)

    merged = None
    if general is not None and gid:
        g = general[general[gid].astype(str).str.strip().isin(ids)].copy()
        if not g.empty:
            merged = g
    if structure is not None and sid:
        s = structure[structure[sid].astype(str).str.strip().isin(ids)].copy()
        if not s.empty:
            if merged is None:
                merged = s
            else:
                merged = merged.merge(s, left_on=gid, right_on=sid, how="outer", suffixes=("_general", "_structure"))

    rows = []
    for _, h in hits.drop_duplicates(subset=[npid_col]).iterrows():
        npid = str(h[npid_col]).strip()
        base = {
            "Plant": plant_name,
            "NPASS_ID": npid,
            "NPASS_source_record": clean_text(" | ".join(str(v) for v in h.tolist() if pd.notna(v))),
            "Source": "NPASS",
        }
        if merged is not None:
            row = merged[merged.astype(str).apply(lambda col: col.str.strip().eq(npid)).any(axis=1)]
            if not row.empty:
                rr = row.iloc[0].to_dict()
                for k, v in rr.items():
                    base[f"NPASS_{k}"] = v
        rows.append(base)

    df = pd.DataFrame(rows)

    # Normalize a few useful fields.
    def first_matching(prefixes):
        for c in df.columns:
            nk = norm_key(c)
            if any(p in nk for p in prefixes):
                return c
        return None

    namec = first_matching(["commonname", "compoundname", "naturalproductname"])
    smic = first_matching(["canonicalsmiles", "smiles"])
    mwc = first_matching(["molecularweight", "mw"])

    if namec:
        df["NPASS_name"] = df[namec]
    else:
        df["NPASS_name"] = ""

    if smic:
        df["NPASS_SMILES"] = df[smic]
    else:
        df["NPASS_SMILES"] = ""

    if mwc:
        df["NPASS_Molecular_Weight"] = df[mwc]
    else:
        df["NPASS_Molecular_Weight"] = ""

    if "NPASS_SMILES" in df.columns:
        rd = df["NPASS_SMILES"].fillna("").map(compute_rdkit)
        rd_df = pd.DataFrame(list(rd))
        rd_df = rd_df.add_prefix("NPASS_")
        df = pd.concat([df.reset_index(drop=True), rd_df.reset_index(drop=True)], axis=1)

    log(f"NPASS: found {len(df)} unique natural products linked to {plant_name}")
    return df


def merge_imppat_npass(imppat_df, npass_df):
    if imppat_df is None:
        return npass_df
    if npass_df is None or npass_df.empty:
        return imppat_df

    # Match first by PubChem CID where available; otherwise by canonicalized SMILES.
    out = imppat_df.copy()
    np = npass_df.copy()

    np_cid_col = None
    for c in np.columns:
        if "pubchem" in norm_key(c) and "cid" in norm_key(c):
            np_cid_col = c
            break

    if "PubChem_CID" in out.columns and np_cid_col:
        out["_join_cid"] = out["PubChem_CID"].fillna("").astype(str).str.replace(r"\.0$", "", regex=True)
        np["_join_cid"] = np[np_cid_col].fillna("").astype(str).str.replace(r"\.0$", "", regex=True)
        merged = out.merge(np.drop_duplicates("_join_cid"), on="_join_cid", how="left", suffixes=("", "_NPASS"))
        return merged.drop(columns=["_join_cid"], errors="ignore")

    # Name fallback is deliberately secondary.
    if "Phytochemical_name" in out.columns and "NPASS_name" in np.columns:
        out["_join_name"] = out["Phytochemical_name"].fillna("").str.lower().str.strip()
        np["_join_name"] = np["NPASS_name"].fillna("").str.lower().str.strip()
        merged = out.merge(np.drop_duplicates("_join_name"), on="_join_name", how="left", suffixes=("", "_NPASS"))
        return merged.drop(columns=["_join_name"], errors="ignore")

    return out


def export_excel(imppat_df, npass_df, merged_df, plant_name, output_path):
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        if imppat_df is not None:
            imppat_df.to_excel(writer, index=False, sheet_name="IMPPAT")
        if npass_df is not None:
            npass_df.to_excel(writer, index=False, sheet_name="NPASS")
        if merged_df is not None:
            merged_df.to_excel(writer, index=False, sheet_name="Combined")
        summary = pd.DataFrame({
            "Field": [
                "Plant",
                "IMPPAT association rows",
                "IMPPAT unique phytochemicals",
                "NPASS matched rows",
                "Application version",
            ],
            "Value": [
                plant_name,
                len(imppat_df) if imppat_df is not None else 0,
                imppat_df["IMPPAT_ID"].nunique() if imppat_df is not None and "IMPPAT_ID" in imppat_df.columns else 0,
                len(npass_df) if npass_df is not None else 0,
                "Version 1.0",
            ],
        })
        summary.to_excel(writer, index=False, sheet_name="Summary")

        for ws in writer.book.worksheets:
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for col_cells in ws.columns:
                max_len = min(max(len(str(c.value or "")) for c in col_cells) + 2, 55)
                ws.column_dimensions[col_cells[0].column_letter].width = max_len


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("850x620")
        self.minsize(760, 540)
        self.log_queue = queue.Queue()
        self._build_ui()
        self.after(150, self._drain_log)

    def _build_ui(self):
        pad = {"padx": 12, "pady": 8}
        frm = ttk.Frame(self)
        frm.pack(fill="both", expand=True, **pad)

        ttk.Label(frm, text="Scientific name of plant:", font=("Segoe UI", 11, "bold")).pack(anchor="w")
        self.plant_var = tk.StringVar(value="Rheum palmatum")
        ttk.Entry(frm, textvariable=self.plant_var, font=("Segoe UI", 12)).pack(fill="x", pady=(0, 10))

        options = ttk.Frame(frm)
        options.pack(fill="x")
        self.imppat_var = tk.BooleanVar(value=True)
        self.npass_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options, text="IMPPAT", variable=self.imppat_var).pack(side="left")
        ttk.Checkbutton(options, text="NPASS", variable=self.npass_var).pack(side="left", padx=15)

        outfrm = ttk.Frame(frm)
        outfrm.pack(fill="x", pady=10)
        ttk.Label(outfrm, text="Excel output:").pack(side="left")
        self.out_var = tk.StringVar()
        ttk.Entry(outfrm, textvariable=self.out_var).pack(side="left", fill="x", expand=True, padx=8)
        ttk.Button(outfrm, text="Browse", command=self.choose_output).pack(side="right")

        ttk.Label(frm, text="Status / log:", font=("Segoe UI", 10, "bold")).pack(anchor="w")
        self.log_box = tk.Text(frm, height=22, wrap="word", state="disabled")
        self.log_box.pack(fill="both", expand=True)

        bottom = ttk.Frame(frm)
        bottom.pack(fill="x", pady=10)
        self.progress = ttk.Progressbar(bottom, mode="indeterminate")
        self.progress.pack(side="left", fill="x", expand=True)
        self.run_btn = ttk.Button(bottom, text="Retrieve & Export Excel", command=self.start)
        self.run_btn.pack(side="right", padx=(10, 0))

        note = (
            "Version 1 retrieves records from the live public database pages. "
            "It does not invent missing phytoconstituents. RDKit descriptors are calculated "
            "from the retrieved SMILES; database-supplied values are kept separately."
        )
        ttk.Label(frm, text=note, wraplength=800, foreground="#555").pack(anchor="w")

    def choose_output(self):
        plant = safe_filename(self.plant_var.get().strip() or "plant")
        p = filedialog.asksaveasfilename(
            title="Save Excel file",
            defaultextension=".xlsx",
            initialfile=f"{plant}_phytoconstituents.xlsx",
            filetypes=[("Excel workbook", "*.xlsx")],
        )
        if p:
            self.out_var.set(p)

    def log(self, msg):
        self.log_queue.put(str(msg))

    def _drain_log(self):
        try:
            while True:
                msg = self.log_queue.get_nowait()
                self.log_box.configure(state="normal")
                self.log_box.insert("end", msg + "\n")
                self.log_box.see("end")
                self.log_box.configure(state="disabled")
        except queue.Empty:
            pass
        self.after(150, self._drain_log)

    def start(self):
        plant = self.plant_var.get().strip()
        if not plant:
            messagebox.showerror("Missing plant", "Enter a scientific plant name.")
            return
        if not self.imppat_var.get() and not self.npass_var.get():
            messagebox.showerror("No source selected", "Select IMPPAT and/or NPASS.")
            return

        output = self.out_var.get().strip()
        if not output:
            output = str(Path.cwd() / f"{safe_filename(plant)}_phytoconstituents.xlsx")
            self.out_var.set(output)

        self.run_btn.configure(state="disabled")
        self.progress.start(12)
        self.log_box.configure(state="normal")
        self.log_box.delete("1.0", "end")
        self.log_box.configure(state="disabled")

        threading.Thread(target=self.worker, args=(plant, output), daemon=True).start()

    def worker(self, plant, output):
        try:
            self.log(f"Starting Version 1 for: {plant}")
            imppat_df = None
            npass_df = None

            if self.imppat_var.get():
                imppat_df = extract_imppat_plant(plant, log=self.log)
                self.log(f"IMPPAT unique phytochemicals: {imppat_df['IMPPAT_ID'].nunique() if not imppat_df.empty else 0}")

            if self.npass_var.get():
                cache = Path(__file__).resolve().parent / "NPASS_cache"
                try:
                    npass_df = extract_npass_from_cache(plant, cache, log=self.log)
                except Exception as e:
                    self.log(f"NPASS automatic retrieval did not complete: {e}")
                    self.log("IMPPAT export will still be produced. You can place official NPASS download files in NPASS_cache and rerun.")
                    npass_df = pd.DataFrame()

            combined = merge_imppat_npass(imppat_df, npass_df)
            export_excel(imppat_df, npass_df, combined, plant, output)

            self.log(f"Excel created: {output}")
            self.log("Finished.")
            self.after(0, lambda: messagebox.showinfo("Completed", f"Excel file created:\n{output}"))
        except Exception as e:
            self.log("ERROR: " + str(e))
            self.log(traceback.format_exc())
            self.after(0, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.after(0, lambda: self.progress.stop())
            self.after(0, lambda: self.run_btn.configure(state="normal"))


if __name__ == "__main__":
    App().mainloop()
