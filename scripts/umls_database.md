# UMLS SQLite Database Setup

This document explains the steps taken to build the **UMLS SQLite database** using `sqlite-utils`.  
These commands are **run in the terminal**, not in Python or R.

---

## 1. Create the SQLite Database

Run once:

```bash
sqlite-utils create-database aux_data/umls.db
```

---

## 2. Import `MRCONSO.RRF`

The `--no-headers` flag is essential because RRF files do not contain header rows.

```bash
sqlite-utils insert aux_data/umls.db MRCONSO raw_data/UMLS/META/MRCONSO.RRF \
    --tsv --delimiter="|" --no-headers --batch-size=50000
```

---

## 3. Rename Key Columns

`MRCONSO.RRF` defines the following key fields:

- Column 1 → **CUI**
- Column 12 → **SAB**
- Column 13 → **TTY**
- Column 15 → **STR**

Rename them:

```bash
sqlite3 aux_data/umls.db "
ALTER TABLE MRCONSO RENAME COLUMN untitled_1 TO CUI;
ALTER TABLE MRCONSO RENAME COLUMN untitled_12 TO SAB;
ALTER TABLE MRCONSO RENAME COLUMN untitled_13 TO TTY;
ALTER TABLE MRCONSO RENAME COLUMN untitled_14 TO CODE;
ALTER TABLE MRCONSO RENAME COLUMN untitled_15 TO STR;
"
```

---

## 4. Create Indexes for Fast Lookup

```bash
sqlite-utils create-index aux_data/umls.db MRCONSO CUI
sqlite-utils create-index aux_data/umls.db MRCONSO SAB
sqlite-utils create-index aux_data/umls.db MRCONSO TTY
sqlite-utils create-index aux_data/umls.db MRCONSO CODE
sqlite-utils create-index aux_data/umls.db MRCONSO STR
```

---

## 5. Test Your Database

Retrieve the first five rows:

```bash
sqlite3 aux_data/umls.db "SELECT CUI, STR, SAB, TTY FROM MRCONSO LIMIT 5;"
```

---

## 6. Import MRREL.RRF

The `MRREL.RRF` file contains relationship information between concepts. Like other RRF files, it does not have headers, so we use the `--no-headers` flag.

```bash
sqlite-utils insert aux_data/umls.db MRREL raw_data/UMLS/META/MRREL.RRF \
    --tsv --delimiter="|" --no-headers --batch-size=50000
```

---

## 7. Rename Key Columns in MRREL

Key columns in `MRREL.RRF` include:

- Column 1 → **CUI1**
- Column 4 → **REL**
- Column 5 → **CUI2**
- Column 7 → **RELA**

Rename them:

```bash
sqlite3 aux_data/umls.db "
ALTER TABLE MRREL RENAME COLUMN untitled_1 TO CUI1;
ALTER TABLE MRREL RENAME COLUMN untitled_4 TO REL;
ALTER TABLE MRREL RENAME COLUMN untitled_5 TO CUI2;
ALTER TABLE MRREL RENAME COLUMN untitled_7 TO RELA;
ALTER TABLE MRREL RENAME COLUMN untitled_11 TO SAB;
"
```

---

## 8. Create Indexes for MRREL

Create indexes to speed up lookups on important columns:

```bash
sqlite-utils create-index aux_data/umls.db MRREL CUI1
sqlite-utils create-index aux_data/umls.db MRREL REL
sqlite-utils create-index aux_data/umls.db MRREL CUI2
sqlite-utils create-index aux_data/umls.db MRREL RELA
sqlite-utils create-index aux_data/umls.db MRREL SAB
```

---

## 9. Import MRSTY.RRF

The `MRSTY.RRF` file contains semantic type information for concepts. Import it similarly:

```bash
sqlite-utils insert aux_data/umls.db MRSTY raw_data/UMLS/META/MRSTY.RRF \
    --tsv --delimiter="|" --no-headers --batch-size=50000
```

---

## 10. Rename Key Columns in MRSTY

Key columns in `MRSTY.RRF` include:

- Column 1 → **CUI**
- Column 2 → **TUI**
- Column 3 → **STN**
- Column 4 → **STY**

Rename them:

```bash
sqlite3 aux_data/umls.db "
ALTER TABLE MRSTY RENAME COLUMN untitled_1 TO CUI;
ALTER TABLE MRSTY RENAME COLUMN untitled_2 TO TUI;
ALTER TABLE MRSTY RENAME COLUMN untitled_3 TO STN;
ALTER TABLE MRSTY RENAME COLUMN untitled_4 TO STY;
"
```

---

## 11. Create Indexes for MRSTY

Create indexes for efficient querying:

```bash
sqlite-utils create-index aux_data/umls.db MRSTY CUI
sqlite-utils create-index aux_data/umls.db MRSTY TUI
sqlite-utils create-index aux_data/umls.db MRSTY STY
```

---

## 12. Import MRSAT.RRF

The `MRSAT.RRF` file contains attribute information for concepts. Import it as well:

```bash
sqlite-utils insert aux_data/umls.db MRSAT raw_data/UMLS/META/MRSAT.RRF \
    --tsv --delimiter="|" --no-headers --batch-size=50000
```

---

## 13. Rename Key Columns in MRSAT

Key columns in `MRSAT.RRF` include:

- Column 1 → **CUI**
- Column 2 → **LUI**
- Column 3 → **SUI**
- Column 4 → **AUI**
- Column 5 → **SATUI**
- Column 6 → **ATUI**
- Column 7 → **CVF**
- Column 8 → **ATN**
- Column 9 → **SAB**
- Column 10 → **ATV**
- Column 11 → **SUPPRESS**
- Column 12 → **CVF2**

Rename them:

```bash
sqlite3 aux_data/umls.db "
ALTER TABLE MRSAT RENAME COLUMN untitled_1 TO CUI;
ALTER TABLE MRSAT RENAME COLUMN untitled_2 TO LUI;
ALTER TABLE MRSAT RENAME COLUMN untitled_3 TO SUI;
ALTER TABLE MRSAT RENAME COLUMN untitled_4 TO AUI;
ALTER TABLE MRSAT RENAME COLUMN untitled_5 TO SATUI;
ALTER TABLE MRSAT RENAME COLUMN untitled_6 TO ATUI;
ALTER TABLE MRSAT RENAME COLUMN untitled_7 TO CVF;
ALTER TABLE MRSAT RENAME COLUMN untitled_8 TO ATN;
ALTER TABLE MRSAT RENAME COLUMN untitled_9 TO SAB;
ALTER TABLE MRSAT RENAME COLUMN untitled_10 TO ATV;
ALTER TABLE MRSAT RENAME COLUMN untitled_11 TO SUPPRESS;
ALTER TABLE MRSAT RENAME COLUMN untitled_12 TO CVF2;
"
```

---

## 14. Create Indexes for MRSAT

Create indexes on key columns to improve query performance:

```bash
sqlite-utils create-index aux_data/umls.db MRSAT CUI
sqlite-utils create-index aux_data/umls.db MRSAT LUI
sqlite-utils create-index aux_data/umls.db MRSAT SUI
sqlite-utils create-index aux_data/umls.db MRSAT AUI
sqlite-utils create-index aux_data/umls.db MRSAT SAB
```

## 15. Overview Diagram — How UMLS Tables Connect

                   +----------------+
                   |   MRCONSO      |
                   |----------------|
 Raw Term  ----->  | STR            |
                   | CUI  <---------------------------+
                   | SAB (MDR/SNOMED/ICD/etc.)        |
                   | TTY (PT/LLT/SY/etc.)             |
                   +----------------------------------+
                                   |
                                   | CUI
                                   |
                                   v
                   +-----------------------------+
                   |          MRREL              |
                   |-----------------------------|
                   | CUI1  --- child concept     |
                   | REL / RELA (ISA, has_SOS)   |
                   | CUI2  --- parent concept    |
                   | SAB = 'MDR' (MedDRA only)   |
                   +-----------------------------+
                                   |
                                   | parent CUI
                                   v
                   +-----------------------------+
                   |          MRSAT              |
                   |-----------------------------|
                   | CUI                         |
                   | ATN (attribute name)        |
                   | ATV (attribute value)       |
                   | → e.g. "System Organ Class" |
                   +-----------------------------+
                                   |
                                   | CUI
                                   v
                   +-----------------------------+
                   |          MRSTY              |
                   |-----------------------------|
                   | CUI                         |
                   | TUI (semantic type)         |
                   | STY (detailed category)     |
                   +-----------------------------+


## 16. Sanity Checks

### 16.1 Test Parent–Child Relationships in MRREL

```bash
sqlite3 aux_data/umls.db "
SELECT c.STR AS term,
       a.ATV AS soc
FROM MRCONSO c
JOIN MRSAT a ON c.CUI = a.CUI
WHERE lower(c.STR) LIKE '%breast cancer%'
  AND a.ATN = 'SOC'
LIMIT 20;
"
```

**Expected:** Parent terms such as *Neoplasms*, *Breast Neoplasms*, *Malignant Neoplasms*.

### 16.2 Test SOC Lookup Through MRSAT

```bash
SELECT c.STR,
       a.ATN,
       a.ATV
FROM MRCONSO c
JOIN MRSAT a ON c.CUI = a.CUI
WHERE a.ATN LIKE '%SOC%'
  AND lower(c.STR) LIKE '%diabetes%'
LIMIT 10;
```

**Expected:**  
`ATV = Metabolism and nutrition disorders`

### 16.3 Test Semantic Type (MRSTY)

```bash
SELECT c.STR, s.STY
FROM MRCONSO c
JOIN MRSTY s ON c.CUI = s.CUI
WHERE lower(c.STR) LIKE '%breast cancer%'
LIMIT 10;
```

**Expected:**  
Semantic types such as *Neoplastic Process* or *Disease or Syndrome*.

