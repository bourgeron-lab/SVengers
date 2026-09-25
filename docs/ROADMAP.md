# Vedette roadmap

État au 2026-09-08 : cœur v0.1 fonctionnel bout-en-bout (call → filter →
merge → genotype → cohort-filter → ped-annotate), validé sur le trio `test`
chr20 (breakpoints identiques à ClinSV, MVR 0,6 %), 34 tests unitaires.

## Prochaines étapes (par priorité)

### 0. Production : models/v5 + seuil SQ>=5 (2026-09-18)
v5 = GIAB + 16 trios wavy GHFC (labels par transmission, `vedette label-trio`).
Holdout 8 enfants wavy : DUP PR-AUC 0,908 (v4) -> 0,988 ; DEL 0,989 -> 0,997.
HG002 : non-régression (SQ5 : P 0,951 / R 0,334 vs v4 0,953 / 0,337).
Filtre optimisé FAIT : scoring batché (2026-09-18, 6h05 -> 21 min) puis
extraction de features avec cache depth par contig (2026-09-21, 21 min ->
13 s/génome avec v5 ; VCF et features TSV identiques à la référence).

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
  les gt.tsv par échantillon, chaque fichier lu une seule fois) ; (j) rapport :
  assembler les matrices de génotypes depuis genotypes/*.gt.tsv plutôt que de
  re-décoder le VCF cohorte (1,4 Md de génotypes ≈ 1 h pysam vs ~5 min TSV).

### 1quater. Garde-fou de génotypage par échantillon (2026-09-25)

Suite du 1ter : le recalcul du garde-fou de novo en excluant `HighCarrierFrac`/
`CNDiscordant` ne changeait rien (déjà mutuellement exclusifs de PASS — 0
site cumule les deux). La vraie cause des plus gros outliers (SMP-05/U/X,
12k-14k candidats chacun) : trois trios **non apparentés** avec un
génotypage globalement peu fiable (60-69&nbsp;% des porteurs PASS à GQ&lt;20,
médiane 10-13) — un problème d'échantillon, pas de site.

Deux garde-fous ajoutés :
- `qc.WAVE_WARN` abaissé de 0,145 à 0,128 (Youden) — capture SMP-06/X (0,139/
  0,140), sous l'ancien seuil.
- Nouveau : fraction de porteurs PASS à GQ&lt;20 par échantillon
  (`report.low_gq_carrier_fraction`, seuil initial 0,40, calibré sur les 3 cas
  connus). Colonne « % porteurs GQ&lt;20 » + KPI dans le rapport.

**Limite constatée à chaud** : sur le run famille FAM-04 (`--known-sites`
pointé sur dataset_paper), les 9 membres sont TOUS flagués (0,59-0,67) malgré
des scores wave normaux (0,10-0,13) et une profondeur normale — `CommonRef`
retire les sites communs faciles du PASS local et ne laisse que des appels
rares/durs, ce qui gonfle mécaniquement la métrique pour tout le monde. Le
seuil 0,40 n'est donc calibré que pour la grande cohorte de référence (sans
`--known-sites`) ; à valider sur le prochain run complet dataset_paper avec
ce garde-fou avant de s'y fier sur des runs familiaux.

### 1ter. v6 — nettoyer les labels par transmission (constat du 2026-09-24)

Validation IGV d'Alexandre sur la famille FAM-04 (v5, sites PASS) :
10 faux / 3 vrais sur les 13 premiers. Les 10 faux sont portés par 8-9
membres sur 9 et sont tous `CohortArtifact` dans dataset_paper (46-100 % de
porteurs sur 2670) ; les DUP fausses ont CN=1 chez 98 % des porteurs.

Conséquence pour l'entraînement : un artefact partagé parent→enfant est
étiqueté **TP par transmission**. Test sur les 16 trios wavy : la règle « DUP
avec dhffc<0,7 » retire 22,6 % des "TP" et 26 % des FP — indiscernables, alors
que sur GIAB elle ne touche presque rien. Les labels wavy v5 sont donc
contaminés par les artefacts familiaux communs, et v5 a appris à les accepter.

À faire pour v6 :
- `vedette label-trio --known-sites <cohort.vcf.gz>` : exclure des labels (ou
  passer en FP) les appels dont le site de référence a REF_CARRIER_FRAC ≥ 0,5,
  et ceux flagués HighCarrierFrac / CNDiscordant dans le run ;
- réentraîner, comparer sur HG002 (non-régression) et sur les validations IGV
  (TSV `validation` exportés des pages events — les croiser par site_id ou
  par locus) ;
- côté post-traitement (fait, 2026-09-24) : `cohort-filter` pose
  HighCarrierFrac / CNDiscordant / CommonRef (`--known-sites`) et le merge
  propage SQMAX/SQMEAN ; le run famille doit pointer `--known_sites` sur la
  dernière grande cohorte.
- nommage : `CohortArtifact` (HWE) attrape aussi de vrais CNV communs
  (site_199 DEL 10,5 kb à 32 % de porteurs, validé vrai) — envisager
  `CommonCNV` distinct pour les sites en déséquilibre HWE mais CN-cohérents.

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
