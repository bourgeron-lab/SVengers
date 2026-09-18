# Vedette roadmap

État au 2026-09-08 : cœur v0.1 fonctionnel bout-en-bout (call → filter →
merge → genotype → cohort-filter → ped-annotate), validé sur le trio `test`
chr20 (breakpoints identiques à ClinSV, MVR 0,6 %), 34 tests unitaires.

## Prochaines étapes (par priorité)

### 0. Production : models/v5 + seuil SQ>=5 (2026-09-18)
v5 = GIAB + 16 trios wavy GHFC (labels par transmission, `vedette label-trio`).
Holdout 8 enfants wavy : DUP PR-AUC 0,908 (v4) -> 0,988 ; DEL 0,989 -> 0,997.
HG002 : non-régression (SQ5 : P 0,951 / R 0,334 vs v4 0,953 / 0,337).
Scoring du filtre batché FAIT (2026-09-18) : 6h05 -> 21 min/génome avec v5
(57 212 records byte-identiques) ; le coût restant est l'extraction de
features par call (accès store), candidate au même traitement par contig.

### 0bis. (archive) Production v4 (2026-09-17)
Pipelines GHFC basculés sur models/v4 (30 features, wave_score),
low_sq_threshold=5 dans configs/ghfc.toml (calibration décalée : SQ5 v4 =
point d'opération SQ7 v3 ; HG002 v5.0q : P 0,953 / R 0,337 / F1 0,498).
Prochain gain attendu : réentraînement enrichi en génomes wavy GHFC.

### 1. Entraînement du filtre ML — v1 FAIT (2026-09-09)
Modèles `models/v1` entraînés sur HG002 40x NovaSeq PCR-free aligné façon GHFC,
vérité GIAB v5.0q stvar, labels Truvari (5644 TP / 60373 FP).
- DEL : PR-AUC(OOF) 0,952. Éval génome vs vérité à SQ≥7 : **précision 99,1 %
  (FDR 0,9 %, cible gnomAD 2,8 % battue ×3), recall 25,9 % global** ; par taille :
  300 pb-1 kb 53 %/99,6 %, 1-10 kb 66 %/99,2 %, 10-100 kb 48 %/100 %.
  Faiblesses : 50-300 pb (13,7 %) et segdup (6,2 % recall à 98,3 % précision).
  Recall brut candidats (pré-filtre) : 40,8 % ≈ Manta.
- DUP : PR-AUC(OOF) 0,643, SQ≥7 → FDR 7,8 % ≈ cible 7,3 %. Vs vérité INS :
  précision 97,2 %, recall 1,4 % (la vérité INS contient surtout des insertions
  non-DUP — l'éval DUP propre exige d'isoler le sous-ensemble tandem-dup de la
  vérité, cf. backlog bench).
- Backlog v1.1 : (a) recall petites DEL 50-300 pb (jonctions SR-only élaguées
  trop fort, matching truvari des DEL en représentation séquence) ; (b) FP DUP
  en amont (nomination everted plus stricte, deux côtés requis) ; (c) HG003/4
  semi-supervisé ; (d) éval DUP tandem-only FAITE (12321/28942 INS tandem ; SQ≥7 : 91,3 % précision vs Manta 60,1 % à recall égal, réservoir brut 5× Manta) ;
  (e) rerun avec PoN fondateurs génome entier ; (f) normaliser la convention
  SVLEN/END des VCF (htslib force END=POS+|SVLEN| sur les allèles symboliques,
  +1 vs notre demi-ouvert — cohérent partout mais à assainir avant publication) ;
  (g) réentraîner v4 avec la feature wave_score ; (h) critère de rescue du
  merge adaptatif à la taille de cohorte (2670 éch. : 80 % des sites sont
  RESCUE avec ΣSQ>=8/2 porteurs — resserrer en fraction de porteurs + ΣSQ
  moyen, ou calibrer sur le taux d'artefacts CohortArtifact) ; (i) workspaces
  familles en une passe batchée (le subset VCF par famille relit tout le VCF
  cohorte : 1207 familles x 635M génotypes = infaisable — assembler depuis
  les gt.tsv par échantillon, chaque fichier lu une seule fois).

### 1bis. Détail du plan d'entraînement initial (référence)
- Télécharger (proxy configuré) dans `resources/giab/` :
  - CRAM HG002/HG003/HG004 ~35x GRCh38 (GIAB FTP, PCR-free ILMN)
  - Vérité HG002-Q100 v1.1 (`.../T2T-Q100_SVs_v1.1/GRCh38/`) + CMRG v1.00 SV
  - GIAB genome-stratifications (segdup/LCR/TR) + UCSC simpleRepeat + Umap k100
- `vedette call` permissif sur HG002/3/4 → labels par `truvari bench`
  (`--pctseq 0 --pctsize 0.7 --refdist 1000`, restreint au BED de confiance ;
  hors-BED = exclu de l'entraînement).
- 2 modèles HistGradientBoosting séparés DEL / DUP, contraintes monotones
  (P(vrai DEL) décroissante en DHFFC, croissante en SR), CV par groupes de
  chromosomes, calibration isotonique out-of-fold, SHAP → élagage features.
- Livraison : `models/v1/` (joblib ×2 + manifest.json {version, schema_hash,
  URLs des données, table seuil→FDR}). Cible : FDR ≤ 2,8 % DEL / 7,3 % DUP
  (= gnomAD-SV v4).
- Augmentation trios GHFC : PAS comme labels v1 (circularité) ; MVR/taux de
  transmission = métriques d'évaluation uniquement.

### 2. Harness de benchmark (`bench/`) — produit les figures du papier
- `bench truth` : Truvari bench+refine vs Q100/CMRG, strates taille
  {50-300, 300-1k, 1-10k, 10-100k, 100k-1M, >1M} × contexte {clean, segdup,
  LCR, TR} ; sortie TSV+JSON.
- `bench compare` : Manta/DELLY/dysgu (installés dans tools/) + WisecondorX/
  ClinSV existants, même protocole ; Vedette avant/après filtre.
- `bench trios` : MVR et transmission en fonction du seuil SQ sur les cohortes
  GHFC (figure de spécificité indépendante de GIAB).

### 3. Hidden split reads — FAIT (2026-09-08)
- `breakpoints/resolve.py` : réalignement fit (numba, sans dépendance) des clips
  stockés (schéma store 0.2), géométrie clip→jonction, micro-homologie HOMLEN,
  dilution sous-bin dans le génotypeur. Validé chr20 : breakpoints ClinSV exacts
  + petites DUP 50-105 pb inédites.
- Reste pour v1.1 : récupération des extrémités riches en mismatches (NM/MD)
  sans clip, et ajout de HOMLEN comme feature du filtre ML (schéma v2).

### 4. Panel of normals — FAIT (2026-09-08)
- `cohort/pon.py` : médiane+MAD par bin en streaming (mémoire plate à 10k
  échantillons), stratifié par sexe, `vedette call --pon`, catalogue
  d'artefacts (`pon artifacts` + `filter --artifact-bed`).
- Sémantique validée sur chr20 GHFC (panel 20 fondateurs EAGER) : segments
  depth 264 → 14, 12/14 survivants rares. Les CNV COMMUNS sortent du canal
  depth (médiane/MAD) et ne restent détectables que par jonctions → --pon
  réservé aux analyses de variants rares, pas aux catalogues d'AF.
- Reste v1.1 : SVD optionnelle (k≈20) façon SavvyCNV ; PoN génome entier
  depuis les ~2000 fondateurs GHFC (array SLURM sur les stores).

### 5. Industrialisation
- Wrapper Nextflow calqué sur tools/ghfc-ngs (DSL2, profils slurm+apptainer).
- INV puis BND/translocations (les évidences PE INV/INTER sont déjà stockées).
- Packaging : PyPI (`vedette` — réservation à faire, ACCORD D'ALEXANDRE requis
  avant toute publication), bioconda, conteneur, mkdocs, CITATION.cff, Zenodo.
- Perf : profiler `extract/scan.py` ; si >40 % du runtime, port Rust ciblé
  (PyO3/rust-htslib) de ce seul module.

## Décisions de design verrouillées (ne pas re-débattre sans raison)
- Merge cohorte = clustering à centroïde borné, PAS de chaînage transitif
  (leçon SURVIVOR : inflation d'AF 1,7-3,6×).
- Re-génotypage cohorte depuis les stores UNIQUEMENT (pattern UKB
  svimmer+GraphTyper2) ; jamais de joint calling complet.
- Filtre = gradient boosting interprétable calibré, PAS de deep learning
  (bioRxiv 2026 Sedlazeck : RF simple = ResNet/Evo2 à F1 95,7 %).
- Vraisemblance de génotype unifiée RD+PE/SR+BAF = la revendication
  méthodologique n°1 du papier ; BAF aussi dans les émissions du HMM.
