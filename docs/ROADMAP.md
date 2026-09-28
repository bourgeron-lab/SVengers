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

### 2bis. Score d'impact plus discriminant (demande Alexandre, 2026-09-28)

Le score actuel est grossier : 0 intergénique / 1 intronique / 2 exonique,
+1 si `candidate_ndd_v5` (7316 gènes), +2 si `hc_ndd_v7` (2200) — soit 4
valeurs possibles et 2 colonnes lues sur les **71** de
`genetrek-data-2026-09-15.tsv`. À ce grain, la moitié du génome code se
retrouve au même niveau.

Colonnes disponibles et directement exploitables (remplissage mesuré sur les
45 045 gènes du fichier) :

- **contrainte / dosage** : `loeuf` (17 961), `shet_mean` (16 486),
  `is_shet_constrained` (3946), `constraint` (4475), `moeuf`, `soeuf`.
  C'est le bon a priori pour une **DEL** (intolérance à la perte de fonction).
  Attention : LOEUF/s_het ne disent rien de la **triplosensibilité** — une DUP
  ne doit pas être notée avec la même échelle (piste : ClinGen triplosensitivity,
  absent du fichier, à ajouter).
- **panels par phénotype** : `panelapp_eng_id_*` (1514), `panelapp_aus_autism_*`,
  `*_epilepsy_*`, `sysndd_definitive_*`, `spark_*`, `sfari_1_*` (244),
  `autism_combined_v7`, `acmg_sf_v3_3` (findings secondaires à signaler).
  Permet une priorisation **selon la cohorte** au lieu d'un score unique.
- **mécanisme** : `chromatin_remodeling`, `epifactors_v2_1`, `syngo_v1_3`,
  `trrust_v2_tf` / `rummagene_tf`, `chrom_epi_tf_combined`.
- **mode de transmission** : `inheritance`, `inheritance_panelapp_*`,
  `inheritance_sysndd_*` — un gène AD touché par une DEL hétérozygote ne pèse
  pas comme un gène AR.
- **évidence quantitative autisme** : `37365347_autism_or`, `_p_value`,
  `_number_of_carriers_among_autistic_individuals`.
- `mane_select` (19 314) pour décider proprement de l'atteinte exonique, et
  `omim_id` pour le lien clinique.

Proposition : score continu séparé par type (DEL/DUP), combinant
(a) la fraction du transcrit MANE touchée, (b) la contrainte du gène le plus
contraint chevauché, (c) l'appartenance au panel du phénotype de la cohorte,
(d) le mode de transmission — plutôt qu'un entier de 0 à 4. Garder un libellé
lisible à côté du nombre, comme aujourd'hui. À caler sur les listes de de novo
déjà validées.

### 2ter. Vérification duphold-like au niveau du génotype (demande Alexandre, 2026-09-28)

Ce que fait duphold (brentp) : pour chaque appel, un fold-change de profondeur
contre les flancs (DHFFC), contre la médiane du génome (DHFC) et contre des bins
de GC comparable (DHBFC), plus le support de jonction.

**On le fait déjà, mais au mauvais endroit.** `filter/features.py` calcule
`dhffc` (exactement DHFFC), `fc_chrom` (~DHFC), `fc_residual`, `depth_cv`,
`frac_consistent`, et le modèle s'en sert — mais **uniquement sur l'appel par
échantillon**. Le re-génotypage de cohorte, lui, ne fait aucune vérification de
ce type : c'est précisément là que les erreurs trouvées le 2026-09-28 se
produisent (96 % des porteurs DUP jugés faux sous IGV ont RD = 1,00, profondeur
strictement normale).

À faire :
- porter le contrôle dans `cohort/genotyper.py` : refuser un génotype porteur
  dont la profondeur contredit le type quand la taille rend la profondeur
  informative (mesuré : de novo apparent 10,1 % pour les DEL à profondeur
  incohérente contre 1,1 % pour les cohérentes, 5,3 % vs 2,1 % pour les DUP) ;
- ajouter le **DHBFC** (comparaison à des bins de GC comparable), le seul des
  trois de duphold qui nous manque — utile là où la correction GC laisse un
  biais résiduel ;
- exposer ces valeurs par génotype dans le VCF cohorte (FORMAT) et dans les
  pages events, pour qu'elles soient lisibles pendant la relecture IGV ;
- corriger en même temps `_junction_support`, qui compte les clips proches du
  point de cassure sans vérifier qu'ils soutiennent *cette* jonction : dans la
  validation, tous les faux porteurs avaient du support de jonction, donc ce
  canal ne discrimine rien en l'état.

Le plancher GQ posé le 2026-09-28 (`--min-carrier-gq`) traite le symptôme le
plus visible ; ces trois points traitent la cause.

### 1octies. Spécificité du génotypeur de cohorte : plancher GQ (2026-09-28)

Mesure à grande échelle sur PMS (77 trios, 32 386 sites échantillonnés) du
**taux de de novo apparent** chez l'enfant porteur quand les deux parents sont
0/0 confiants — un proxy direct du taux d'erreur de génotypage :

| type | sous-groupe | porteurs | de novo apparent |
|---|---|---|---|
| DEL | GQ < 20 | 20 832 | **16,7 %** |
| DEL | GQ >= 20 | 251 294 | **0,8 %** |
| DUP | GQ < 20 | 71 508 | **13,9 %** |
| DUP | GQ >= 20 | 253 823 | **1,9 %** |
| DEL | profondeur incohérente | 27 069 | 10,1 % |
| DUP | profondeur incohérente | 254 989 | 5,3 % |

Le GQ domine tous les autres critères testés (cohérence de profondeur, taille,
support de jonction — ce dernier ne discrimine rien, tous les faux en ont,
signe que le comptage de reads clivés est trop permissif : il compte les clips
proches du point de cassure sans vérifier qu'ils soutiennent *cette* jonction).

Fait : sur les pages events, « porteur » = génotype alternatif **et**
`GQ >= --min-carrier-gq` (défaut 20). Les génotypes moins sûrs restent
affichés dans les colonnes par membre mais ne comptent plus comme portage.

Vérification contre les 42 sites relus sous IGV (famille FAM-04) :

| critère | vrais | faux | précision |
|---|---|---|---|
| brut | 3/3 | 39/39 | 7 % |
| GQ>=20 | 3/3 | 18/39 | 14 % |
| adossé à un appel propre | 3/3 | 18/39 | 14 % |
| **GQ>=20 ET adossé** | **3/3** | **10/39** | **23 %** |

Aucun vrai perdu, 29 faux sur 39 écartés. Reste à traiter :
- la spécificité du comptage SR/PE dans `_junction_support` (apparier les
  clips à la jonction, pas seulement à la position) ;
- les statistiques de cohorte (NCARRIER/AF/CARRIER_FRAC) comptent toujours
  tous les génotypes alternatifs, sans plancher GQ — à trancher, car y
  toucher déplace AF, HWE et les tags dérivés.

### 1septies. Le portage par re-génotypage seul est quasi toujours faux (2026-09-28)

Validation IGV **en aveugle** d'Alexandre sur 60 événements de la famille
FAM-04 dans PMS (40 sites libérés par le critère directionnel + 20
contrôles déjà PASS, mélangés, FILTER masqué), 42 relus :

| groupe | n | vrai | faux |
|---|---|---|---|
| libérés | 25 | 3 | 22 |
| contrôles PASS | 17 | 0 | 17 |

Les contrôles à 0/17 disent que le problème préexiste au changement
directionnel. La cause est ailleurs, et elle est nette — croisement avec les
appels **propres** des membres de la famille :

| verdict | appel propre PASS | appel propre filtré | aucun appel propre |
|---|---|---|---|
| vrai (3) | 3 | 0 | 0 |
| faux (39) | 6 | 12 | 21 |

Un site de cohorte est un groupe de points de cassure ; le re-génotypage
attribue ensuite un génotype à **tous** les échantillons, qu'ils aient produit
un appel ou non. Sur la famille FAM-04 dans PMS : **80,4 % des 19 416 événements
PASS portés n'ont aucun appel propre derrière eux** (11,8 % ont un appel PASS,
7,8 % un appel filtré). La relecture note 0/33 pour les non-adossés contre
3/9 pour les adossés PASS.

Cela explique l'écart avec le benchmark HG002 (précision 0,95 à SQ>=5) : celui-ci
mesure les **appels par échantillon**, pas le portage re-génotypé affiché dans
les pages events. Et cela explique pourquoi la même famille paraissait bonne
dans son run isolé (9 échantillons : peu de sites d'autrui à hériter).

Fait : colonne « appel propre » (PASS / filtré / aucun) et sélecteur sur les
pages events, réglé par défaut sur « avec appel propre » — 19 416 -> 3 806
lignes sur la famille FAM-04. Le défaut retient les appels filtrés, sinon les
sites rescués (la DEL SHANK3 familiale, LowSQ chez chaque porteur) seraient
masqués.

À faire : regarder la spécificité du génotypeur de cohorte lui-même (seuil
d'appel d'un porteur au re-génotypage) — c'est la source du problème, la
colonne ne fait que le rendre visible.

### 1sexies. Discordance de génotypage rendue directionnelle (2026-09-28)

Sur la cohorte PMS (Phelan-McDermid, 319 génomes), la délétion SHANK3 causale
— chr22:50 654 001, 33 kb, celle de la famille FAM-04 — sortait `CohortArtifact`.
Ni la fréquence ni HWE (p=1,0) : c'est `GTR_DISC`. 6 échantillons seulement
l'avaient nominée au calling (SQ faible, rescue), le re-génotypage y trouve
64 porteurs sur 319, écart relatif 9,7 > seuil 0,5.

Or c'est le fonctionnement voulu du pipeline : le re-génotypage depuis les
stores est *fait* pour être plus sensible que le calling par échantillon. Le
critère punissait ce gain. Asymétrie mesurée sur PMS : 36 328 sites filtrés
dans le sens « le re-génotypage en trouve plus » contre 425 dans l'autre.

Critère désormais **directionnel** : on ne flague que
`(nominations - porteurs) / nominations > seuil`, c'est-à-dire des appels que
le re-génotypage ne reproduit pas. Le cas artefact inverse (« le génotypeur
appelle tout le monde ») est couvert depuis le 2026-09-24 par
`HighCarrierFrac`, qui n'existait pas quand `GTR_DISC` a été écrit.

Effet sur PMS : 7 486 sites deviennent PASS (+16 %). Comparaison Mendel
(conditionnelle au fait que l'enfant porte, cf. ci-dessous) : MVR 22,9 % pour
les libérés contre 22,9 % pour les PASS actuels et 3,0 % pour les autres
artefacts — ils se comportent comme le jeu accepté, pas comme la classe
artefact (dont 92 % sont portés par tout le monde). Réserve : leur taux de
de novo apparent est plus haut (14,5 % contre 9,7 %). Validation IGV en
aveugle en cours pour trancher.

**Piège de lecture de `mendel_report.tsv`** : le fichier ne contient que les
lignes « informatives » (`mendel.py`, « only report informative rows »), donc
les enfants non-porteurs en sont absents. Tout taux de transmission calculé
directement depuis ce TSV vaut ~1,0 par construction — utiliser les colonnes
du résumé ou `bench trios`, qui comptent sur l'ensemble des sites génotypés.

### 1quinquies. `lcr_ovl` était morte depuis le début (2026-09-25)

Alexandre : « les évènements dans un LCR sont parfois plus faux, surtout les
DUP ». Vérification : `lcr_ovl > 0` pour **0 appel sur 38 094** (génome
entier) et **0 sur 29 056** (HG002, données d'entraînement). Cause : le track
`selfchain-LCR.hg38.bed` est nommé à l'Ensembl (`1`, `2`…) alors que les
appels sont en UCSC (`chr1`) — aucun chevauchement ne pouvait matcher.
`IntervalSet` résout désormais les deux conventions (test dans
test_reference.py).

Conséquences :
- **Aucun impact sur v5** : la feature étant constante à 0 à l'entraînement,
  aucun arbre ne splitte dessus — vérifié, SQ strictement identiques sur
  5000 appels avec `lcr_ovl` aléatoire (écart max 0,0).
- **Pour v6** : la feature devient réellement informative ; à réentraîner
  pour qu'elle serve. C'est aussi une leçon de contrat : une feature
  constante devrait échouer bruyamment — **garde-fou ajouté** le 2026-09-25 :
  `train.check_feature_variance` refuse d'entraîner sur une feature de
  variance nulle (ou entièrement NaN), avec `--allow-constant-features` pour
  forcer. Exemptions par type : `is_dup` (un modèle par type) et, sur DUP,
  les canaux propres aux DEL (`cig`, `ins_z`, `ins_pairs`). Vérifié sur les
  données v5 : le garde-fou aurait refusé l'entraînement en nommant
  `lcr_ovl`. Les features constantes détectées sont aussi tracées dans le
  manifeste du modèle.

Hypothèse LCR × DUP elle-même : **non confirmée** sur les 45 validations
IGV disponibles une fois la feature réparée — DUP faux `lcr>0` 3/8
(médiane 0,000) contre DUP vrais 5/10 (médiane 0,008). Ce qui sépare
vraiment, sur ce petit échantillon : le type (DEL 26 vrai/1 faux = 96 % vs
DUP 10/8 = 56 %) et `segdup_ovl` (3/8 des DUP faux à 1,00, 0/10 des vrais).

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
- Nouveau (corrigé deux fois le 2026-09-25) : **charge de porteurs**
  `report.carrier_excess` = nb de sites **PASS** portés / médiane cohorte,
  seuil **5x**, appliqué seulement à partir de 20 échantillons (critère
  relatif). Distribution dataset_paper : médiane 1,00, p90 4,4x, p95 6,4x,
  p99 9,6x, max 23,5x. Le seuil 5x flague 212/2670 (7,9 %) et capture 18 des
  20 pires outliers de novo (charge médiane des outliers 1,41x contre 0,91x
  pour les autres).

  **Deuxième erreur de calibration** (même classe que la première) : le seuil
  1,3x avait été calibré sur le nombre de porteurs **tous sites** (médiane
  308 177) alors que le garde-fou mesure les porteurs **PASS** (médiane
  3 284) — deux distributions sans rapport, d'où 1037/2670 flagués au premier
  essai. Règle à tenir : calibrer sur *exactement* la quantité que le code
  calcule, en lisant la colonne produite par le code, jamais une quantité
  voisine.

  **Première version erronée** : la fraction de porteurs PASS à GQ&lt;20,
  seuil 0,40 — calibration à l'échelle réelle : 2670/2670 flagués. La
  médiane cohorte de cette fraction est 0,75 (p5 0,68, max 0,89) et les trois
  vrais coupables (SMP-05/U/X) sont à 0,61-0,69, soit les percentiles 2-6 :
  *parmi les plus propres*. Ce qui les sépare, c'est la charge (487-497k
  porteurs contre 308k en médiane, percentile 100). Leçon : ne pas calibrer
  un seuil sur 3 cas sans la distribution de fond. La fraction GQ reste
  affichée (elle décrit le callset : 75 % des génotypes porteurs PASS sont
  sous GQ 20) mais ne sert plus de critère.

Note : l'hypothèse « c'est `CommonRef` qui gonflait la métrique sur la
famille FAM-04 » était fausse — la fraction est haute partout, y compris sur la
cohorte de référence sans `--known-sites`. Le critère relatif actuel reste en
revanche inapplicable aux petits runs familiaux (tout le monde porte les mêmes
sites), d'où le minimum de 20 échantillons.

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
