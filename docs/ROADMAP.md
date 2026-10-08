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

### 1undecies. Bornes des événements familiaux : celles du porteur, pas du représentant (2026-09-30)

Alexandre, sur FAM-10 (délétion SHANK3) : IGV place la délétion à
chr22:45,145,667-50,808,468, la page familiale affichait 44,094,543-50,800,002.

Le caller par échantillon avait juste : l'appel propre du proposant commence à
**45,146,001**, à 300 pb d'IGV. Deux erreurs sont venues ensuite, au report :

1. **Coordonnées du représentant de cluster.** Au merge, l'appel du proposant
   (45,15-48,99 Mb) est entré dans un site dont le représentant — un autre
   patient PMS, la cohorte réunit des 22q13 de toutes tailles — commence à
   44,116,001. Le re-génotypage lui donne 0/1 GQ 99 (79 % du site est bien
   délété chez lui), et la ligne héritait des bornes de l'autre patient : 1 Mb
   de trop. Même défaut le même jour sur **FAM-09** (IGV 43,784,666,
   page 43,570,001), mais là l'appel propre revient en neuf morceaux dont
   aucun ne recouvre le site à 50 % : une règle « un appel propre » ne suffit pas.

   Règle finale (`own_span`) : les appels propres des porteurs qui recouvrent
   le site (≥ 10 % de leur longueur dans le site) sont chaînés avec la règle
   du collapse ; la meilleure chaîne remplace les bornes du site si elle le
   recouvre à RO ≥ 0,5. Les bornes de la chaîne viennent uniquement d'appels
   **PASS, non emboîtés dans un appel plus long, et dont le site de cohorte
   est PASS** ; les autres appels propres servent seulement de liaison.
   Chaque garde répond à un cas PMS réel :
   - appels filtrés comme bornes → 57 656 lignes déplacées sur la foi
     d'appels rejetés (d'où « PASS seulement ») ;
   - FAM-08 : le milieu (138 kb) de la délétion est un appel LowSQ ;
     sans liaison par les appels filtrés, l'événement se coupait en deux ;
   - FAM-11 : la duplication est entrée dans un site commençant 141 kb
     à l'intérieur ; une règle « jamais au-delà du site » la coupait en deux
     (la chaîne peut donc agrandir la ligne) ;
   - FAM-14 : un appel PASS subtélomérique de 23 kb, dont le site est
     porté par les 319 échantillons, poussait la fin au-delà du chr22 (d'où
     « site de cohorte PASS ») ; un appel qui ne touche le site que sur 1 pb
     faisait de même (d'où les 10 %) ;
   - FAM-06 : cinq appels PASS de 0,1-0,8 kb emboîtés dans une
     duplication LowSQ de 712 kb devenaient des bornes et coupaient une
     duplication de 5,6 Mb (d'où « non emboîtés »).

   Résultat sur les 120 familles PMS, contre le rapport du 30/09 12 h 24 :
   FAM-09 passe à 43,785,001 (335 pb d'IGV), FAM-10 / FAM-11 / FAM-15 /
   FAM-01 / FAM-08 / FAM-14 / FAM-06 inchangées, aucune ligne au-delà de
   la fin d'un chromosome ; 3 % des lignes bougent (médiane 1,6 % de leur
   taille) ; parmi les événements ≥ 1 Mb, six changent de borne, chacun vers
   l'appel propre PASS du porteur (FAM-02, FAM-07, FAM-09, FAM-12,
   FAM-13, FAM-05 à 1 pb).
2. **Pont sur un petit voisin.** Un polymorphisme de 0,6 kb (NCAR 38, mère à
   GQ 11 donc non-porteuse) à 21 kb du bord a été collé par la tolérance aux
   trous (2 % de 4,9 Mb = 97 kb). Correctif : franchir un trou exige que le plus
   petit des deux voisins soit au moins aussi long que le trou ; les fragments
   d'un même événement sont grands, les polymorphismes voisins non.

Fin de l'événement : 50,800,002 contre 50,808,468 dans IGV. Ce n'est pas un
écart de bornes : le segment 50,800-50,819 kb est subtélomérique (LowMap, MaxN),
porté par les 319 échantillons, donc non-PASS. 50,808,468 est la fin du
chromosome moins le gap télomérique de 10 kb.

Reste ouvert, au génotypeur cette fois : un porteur partiel reçoit GQ 99 sur un
site qui dépasse son événement de 21 %. Un test de CN sur les bords du site
(premier et dernier 10 %) signalerait ce recouvrement partiel pour les porteurs
**sans appel propre PASS**, que le correctif du report ne peut pas couvrir.

### 1septdecies. `wave_score` : v5 inonde les génomes ondulés ; v5n ; v8 non retenu (2026-10-08)

**Pistes de production reproduites** pour GIAB : `by1000.bedgraph` =
`mosdepth -n -m -x -Q 40 --by 1000` (médiane, mode rapide, MAPQ ≥ 40 :
64 445/64 445 bins identiques sur le chr20 d'un génome GHFC) ; `vaf.bedgraph`
= DeepVariant 1.9.0 (`run_deepvariant --model_type WGS`, paramètres = ceux de
`ghfc-ngs`), sites PASS, QUAL ≥ 30, DP ≥ 10 (0 en trop, 0 manquant).

**Découverte en route : v5 a appris une marche sur `wave_score`** (feature
par génome) entre 0,10 et 0,13. Mêmes appels, mêmes features, seule la wave
change : génome PMS 4 794 → 29 444 appels à SQ ≥ 5 (wave 0,10 → 0,13) ;
HG002 pistes de production (wave 0,164) 42 082 ; un génome dataset_paper
forcé à 0,16 : 5 070 → 42 764. Production touchée :
- PMS (stores schéma 0.4, wave présente) : 102/319 génomes ≥ 0,12, PASS
  médian 12 600 contre 5 500 ;
- dataset_paper : 2 517 stores en schéma 0.3 (wave absente = −1, d'où un v5
  « sain » par accident), mais 153 en 0.4 dont **97 ≥ 0,12** — très
  probablement les ~81 génomes « charge ≥ 5× » des rapports.
HG002 avec pistes de production, v5 : DEL P 0,456 / R 0,327.

**v5n** (`models/v5n`, commit 21ab8bc) = estimateurs v5 inchangés, `wave_score`
neutralisée à −1 au scoring (`neutral_features` du manifeste ; `vedette train
--neutralize`). Identique à v5 sur tout génome sans wave ou sous la marche.
PASS par génome : génomes normaux identiques ; génomes ondulés 10-18k → 3,3-6k
(DUP de jonction < 1 kb 10-11k → 0,4-0,9k). HG002 pistes de production :
DEL P 0,859 / R 0,367 / F1 0,515 (v5 : 0,381) ; CMRG DEL F1 0,563.

**v8 non retenu.** v5-like + GIAB sur pistes de production + labels wavy
nettoyés + `wave_score` neutralisée. Par rapport à v5n :
- PASS par génome 2,8-5,4k (v5n 3,3-5,6k) ; DUP de jonction < 1 kb divisées
  par ~2 ; DEL de profondeur ≥ 1 kb stables (pas l'explosion de v7), sauf un
  génome PMS ondulé (457 contre 119, 17 % récurrentes) ;
- holdout wavy nettoyé : DEL P 0,970 / R 0,585 (v5n 0,829 / 0,635) ; DUP
  P 0,802 / R 0,163 (v5n 0,179 / 0,831) ;
- HG002 hors échantillon, pistes de production : DEL P 0,898 / R 0,290 /
  F1 0,439 contre v5n 0,515 ; CMRG DEL F1 0,472 contre 0,563.
Plus précis, mais le rappel perdu sur HG002 hors échantillon et CMRG
l'emporte. À creuser avant un v9 : la part de la neutralisation de la wave
(qui portait peut-être un vrai signal pour les génomes ondulés) et celle des
labels GIAB refaits.

### 1sexdecies. Second avis par CNN sur les pistes du store — PLAN (2026-10-06)

Question d'Alexandre : faire pour les DEL/DUP ce que DeepVariant fait pour
les SNV/indels, c.-à-d. classer une « image » de l'événement.

**Constat de départ.** Le filtre (gradient boosting, `filter/schema.py`) ne
voit que des résumés : `dhffc`, `depth_cv`, `frac_consistent`,
`n_bins_informative`, `flank_asym`… jamais la *forme* du signal. Cas d'école
(PMS, 2026-10-06) : délétion de 9 kb de l'extrémité 3' de SHANK3 chez un
proposant, CN 1, RD 0,50, GQ 99, unique dans la cohorte, confirmée dans IGV
en une seconde (deux marches nettes à 0,5, flancs plats à 1,0) — SQ 0,38 pour
un seuil de 5 : profondeur seule, aucun site hétérozygote (22q13.33 riche en
GC), `depth_cv` 0,30. À l'inverse une vague de GC a une forme ondulante sans
marche que les résumés distinguent mal (cf. les ~70 vagues CN ~1,83 du scan
mosaïque de dataset_paper, 1quaterdecies). Le palier « sous le seuil, gènes
NDD » des pages familles (2026-10-06) rattrape ce cas pour la relecture ; le
CNN viserait à ce que le score lui-même le voie.

**Littérature.** Samplot-ML (Belyeu et al., GigaScience 2021 : CNN sur images
samplot de DEL, filtrage de faux positifs) ; Cue (Popic et al., Nat Methods
2023 : images multi-canaux des signaux d'alignement, détection d'objets) ;
CNV-espresso (exome : images de profondeur de l'échantillon *et* d'échantillons
de référence, confirmation in silico) ; DeepSV ; SVision (lectures longues).

**Ce qu'on NE fait PAS : rendre des PNG.** Le signal d'un CNV de quelques kb
n'est pas à la base près comme celui d'un SNV : il est déjà condensé dans le
store (`depth`, `baf`, `pe`, `sr`, `ins`, `cig`). Le store *est* l'image ;
rastériser perd de l'information et coûte cher. Invariant respecté : le
second avis ne lit que le store, jamais le CRAM.

**Plan.**
1. *Tenseurs* (`vedette/filter/tensor.py`, nouveau) : pour chaque appel,
   fenêtre = appel + un flanc de même longueur de chaque côté (bornée, min
   10 kb), rééchantillonnée à L = 256 points ; canaux :
   - log2 du ratio de profondeur normalisé diploïde, et masque ;
   - BAF : densité de sites hétérozygotes et médiane |BAF − 0,5| par point ;
   - paires discordantes (DEL/DUP orientées), split reads, décalage de taille
     d'insert (`ins`), indels CIGAR ;
   - **contexte cohorte** : médiane et MAD du ratio de profondeur des génomes
     de la cohorte (ou d'un panel de référence) sur les mêmes bins — le levier
     contre les artefacts systématiques (creux chr16:28-31 Mb de toute la
     cohorte PMS, vagues GC). Coût : un résumé par bin calculé une fois par
     cohorte (comme `vedette mosaic`, ~1-4 min), pas par appel ;
   - canaux scalaires concaténés après la convolution : log10 taille, type,
     source, `popcov_stdev`, `segdup_ovl`, `lcr_ovl`.
   Extraction ~ms par appel ; à borner aux appels pré-filtrés (hard filters
   passés), pas aux ~35k appels bruts par génome.
2. *Modèle* : CNN 1D léger (3-4 blocs conv + pooling, < 1M paramètres), un
   par type ou un commun avec le type en entrée ; CPU suffisant à
   l'inférence. Sortie = `SQ_CNN` calibré (isotonique), écrit en FORMAT à
   côté de `SQ`. Le schéma de features GBM n'est pas touché (append-only).
3. *Labels* — le point dur, cf. v6/v7 (1quindecies) :
   - vérité GIAB HG002/3/4 (Truvari), MAIS avec une profondeur produite comme
     en production (bedgraph 1 kb façon GHFC) — prérequis commun avec v8 ;
     sinon on rejoue le décalage CRAM/bedgraph de v7 ;
   - transmission trio nettoyée (`label-trio --known-sites`) ;
   - **validations IGV des pages familles** (vrai/faux/incertain, export TSV) :
     labels humains sur exactement les cas difficiles. Objectif avant
     entraînement : quelques centaines, dont les appels « sous le seuil ».
   - augmentation : symétrie gauche/droite, léger jitter des bornes, bruit de
     profondeur ; pas d'événements synthétiques injectés comme seule source.
4. *Usage en deux temps* :
   a. **second avis sur les appels LowSQ seuls** (comme SHANK3) : un appel
      rejeté par SQ mais SQ_CNN élevé devient « RESCUE_CNN », au rayon
      d'action limité ;
   b. si (a) tient : combinaison SQ/SQ_CNN (empilement) pour tous les appels.
5. *Critères de déploiement* (leçon v6, mémoire « contrôle avant
   déploiement ») : PASS par génome **par type, source et taille** vs v5 sur
   ≥ 5 génomes dont des wavy ; récurrence des nouveaux PASS chez des
   non-apparentés ; HG002 hors échantillon ; holdout nettoyé ; et, pour (a),
   nombre de rescues par génome (cible ≲ 1) relus dans IGV sur un
   échantillon.

**Estimation** : tenseurs + contexte cohorte ~2 j ; entraînement/évaluation
~2 j une fois les labels IGV disponibles et v8 (profondeur façon production
pour GIAB) faite. Ne pas lancer avant ces deux prérequis.

### 1quindecies. Décalage entraînement/production sur `gq` — à corriger par v6 (2026-10-05)

Benchmark HG002 de la v1.0.0 (docs/benchmark_hg002_v1.tsv) contre celui du
modèle v5 (09-18) : DEL 50-300 pb rappel 0,225 → 0,210, LCR 0,192 → 0,172,
F1 global 0,495 → 0,483 ; CMRG identique. Autopsie des 132 DEL vraies
50-300 pb perdues : 103 sont toujours appelées (JUNCTION) mais tombent en
LowSQ, 21 manquent à l'appel. Le filtre n'y est pour rien (anciens appels
refiltrés avec le code v1 : même SQ médian 5,16). Entre les deux appels, **une
seule feature change pour les 104 candidats appariés : `gq` (+10 en
médiane)** ; remettre l'ancien `gq` sur les features v1 rend PASS 95 des 104
(SQ médian 4,35 → 5,16).

Origine : `fd1f49e` (09-09 17:34) remplace la vraisemblance de Poisson du
nombre de lectures de jonction par une binomiale négative. Les appels GIAB
HG002/3/4 qui ont servi à entraîner v3 → v5 datent du 09-09 matin (Poisson) ;
tous les appels de production (cohortes GHFC, HG002 v1) sont postérieurs
(binomiale négative). Le modèle v5 sert donc en production une feature dont
la distribution a changé depuis son entraînement. Le benchmark v5 (sur les
appels du 09-09) surestimait légèrement la production ; celui de la v1 est le
fidèle.

**v6 entraîné le 2026-10-05** (`models/v6`, pas encore en production) : mêmes
données et labels que v5, appels GIAB refaits avec le code v1, features avec
`lcr_ovl` vivante. Holdout 8 enfants wavy (features de production, hors
entraînement des deux modèles), au seuil SQ ≥ 5 : DEL P 0,974 / R 0,714 (v5)
→ P 0,978 / R 0,973 (v6) ; DUP P 0,795 / R 0,728 → P 0,878 / R 0,960 —
le décalage de `gq` coûtait donc ~25 points de rappel en production. HG002
avec un modèle v6 entraîné SANS HG002 (premier chiffre hors échantillon :
v3-v5 étaient évalués sur un HG002 vu à l'entraînement) : DEL R 0,365 /
P 0,870 / F1 0,514 (v5 en échantillon : 0,323 / 0,955 / 0,483 ; Manta
0,376 / 0,934 / 0,536 ; dysgu 0,389 / 0,879 / 0,539) ; meilleur F1 de la
table dès 300 pb (300 pb-1 kb 0,722, 1-10 kb 0,831, 10-100 kb 0,706) ; CMRG
F1 0,573 (v5 0,531). docs/benchmark_hg002_v6.tsv.

**v6 retiré de la production le 2026-10-06** (mis en production le 05/10,
annulé avant tout rapport). Au seuil SQ ≥ 5, v6 laisse passer 12-17k appels
par génome dataset_paper contre 4,7-5,6k avec v5 — surtout des petites DUP de
jonction < 1 kb (5-8,6k par génome) — et le merge passe de 523k à 1,92M sites
(cohort-filter tué à 16 Go). Le holdout wavy (P 0,878 / R 0,960 en DUP) était
trompeur : il ne note que les appels à label de transmission net et écarte
les ~25k appels « ambigus » par enfant, là où se trouvent ces DUP. HG002 hors
échantillon le disait déjà (246 FP DUP contre 46, précision DEL 0,870).
Leçon : un modèle se valide sur le NOMBRE de PASS par génome et sur les
métriques trios de la cohorte (MVR, de novo apparents) avant la production,
pas seulement sur un holdout aux labels filtrés.
À faire pour un v6 utilisable : (1) nettoyer les labels de transmission
(un artefact récurrent présent chez un parent est « transmis » ; exclure les
sites communs via `--known-sites`, cf. 1ter) ; (2) calibrer le seuil sur une
précision cible mesurée sur l'ensemble des appels (HG002 hors échantillon +
taux de de novo apparents), pas sur les labels nets ; (3) garde-fou au
déploiement : comparer le nombre de PASS par génome au modèle en place.

**v7 (2026-10-06), non retenu non plus.** v6 + labels de transmission
nettoyés (`label-trio --known-sites` : pas de label sur les sites communs ou
marqués artefact ; 99 % des petites DUP de jonction « TP » y tombaient).
Bilan sur quatre critères :
- PASS par génome (6 génomes dataset_paper) : génomes normaux +12 à +26 %
  vs v5 (v6 : ×2,5 à ×3,4) ; génome wavy surchargé −40 % ; autre wavy +66 %.
- holdout wavy (labels nettoyés, sites rares/modérés) : DEL P 0,938 / R 0,870
  (v5 0,829 / 0,635) ; DUP rares faibles pour tous (v5 P 0,18, v7 0,27 à
  rappel divisé par deux).
- HG002 hors échantillon : DEL P 0,883 / R 0,358 / F1 0,509 (v6 0,870 /
  0,365 / 0,514) ; en échantillon v7 0,945 / 0,345 / 0,506 vs v5 0,955 /
  0,323 / 0,483 ; DUP (espace INS) précision 0,939 vs 0,885 ; CMRG F1 0,573.
- MAIS les PASS supplémentaires sur les génomes de cohorte sont des **DEL de
  profondeur seule ≥ 1 kb ×6 à ×18** (269 → 1679 ; 169 → 3142 chez un wavy),
  sans support de jonction, 45 % présentes chez les 3 non-apparentés testés
  (v5 : 28 %) — artefacts de profondeur récurrents très probables. Invisible
  sur HG002 : la profondeur GIAB vient du CRAM (giab.toml), celle des
  cohortes du bedgraph 1 kb précalculé. Une fois les labels wavy (seuls
  exemples « bedgraph ») nettoyés, le modèle apprend la profondeur sur les
  CRAM et la sert sur des bedgraphs : nouveau décalage entraînement/production,
  sur les features de profondeur cette fois.
- Le test de transmission ne tranche pas (90 % présentes chez un parent),
  pour la même raison que les labels : un artefact récurrent est « hérité ».

Pour un modèle utilisable (v8) : produire pour HG002/3/4 la même profondeur
que la production (bedgraph 1 kb façon GHFC, appel en fast path), refaire
features et labels, réentraîner ; critères de déploiement = PASS par génome
**par type et par source** vs v5, récurrence chez des non-apparentés,
HG002 hors échantillon, holdout nettoyé.

Correctif (fait) : **v6** — rappeler HG003/HG004 avec le code v1 (HG002 l'est déjà,
`resources/giab/hg002_v1/`), recalculer les features, mêmes labels, et
réentraîner (avec `lcr_ovl` désormais vivante et des grands événements
réalistes, cf. 1ter / 1nonies) ; re-benchmarker. Garde-fou à ajouter : le
manifeste du modèle doit stocker des quantiles de référence par feature, et
`vedette filter` avertir quand la distribution servie s'en écarte (test de
dérive) — ce décalage serait sorti au premier run de production.

### 1quaterdecies. Mosaïques chromosomiques (`vedette mosaic`, 2026-10-02)

Vérifié par Alexandre dans IGV sur le père SMP-02 (famille FAM-03) : une seule
délétion chr10:83,36 Mb-fin là où Vedette rendait 57 fragments CN=1. La
profondeur donne **CN ≈ 1,33 homogène** sur 50 Mb : une perte dans ~2/3 des
cellules. Le même clone (fraction 0,62-0,67) porte 12 segments ≥ 5 Mb — pertes
2q, 6q, chr9, 10q, 12p, 17p, gains 4p et 17q (profil d'isochromosome 17q) —
avec le déséquilibre allélique attendu (|BAF − 0,5| 0,20 sur 10q, 0,05 chez la
mère ; attendu 0,25 pour une perte, 0,125 pour un gain à f = 0,67). Non
transmis aux enfants : clone somatique (artefact de culture d'une lignée
lymphoblastoïde, ou processus clonal dans le prélèvement), pas un événement
constitutionnel — signalé à Alexandre comme possible découverte fortuite.

Pourquoi Vedette le ratait : le CN-HMM n'a que des états entiers ; un rapport
de 0,665 alterne entre CN 1 et CN 2 (fragments), et les gains à CN 2,67 ne
tombent sur aucun état (non appelés).

Passe ajoutée (`src/vedette/mosaic.py`, `vedette mosaic --run-dir`, ~1 s par
échantillon depuis le store) : fenêtres de 1 Mb, segmentation binaire, fraction
cellulaire sous hypothèse de changement d'une copie (perte CN = 2 − f, gain
CN = 2 + f), puis deux filtres : **homogénéité** (≥ 75 % des fenêtres à ±0,2 du
CN du segment — un mélange de niveaux entiers n'est pas un clone) et **BAF**
(excès de |BAF − 0,5| ≥ 40 % de l'attendu pour f). Mosaïque si
0,15 ≤ f ≤ 0,85 ; les segments à f ≈ 1 (grands gains constitutionnels, ex.
duplications 8q24 des translocations PMS) sont écrits mais sans alerte. Le
rapport lit `samples/<S>/<S>.mosaic.tsv` : alerte « mosaïque 67 % chr10:83-133
Mb » sur les lignes qui tombent dans un segment du porteur, « porteur
mosaïque » sur ses autres lignes, table « Mosaïques chromosomiques ».

Calibrage du filtre d'homogénéité sur PMS : sans lui, une duplication 16p13
constitutionnelle de 3 Mb voisine d'une fenêtre délétée (CN par fenêtre 2,9 /
2,5 / 3,0 / 1,15) et une délétion 22q13 au bord intermédiaire (1,66 / 1,65 /
1,03 / 1,02) passaient pour des mosaïques ; avec lui, seul reste un proposant
à CN 1,60-1,66 homogène sur 22q13 (43-50 Mb) avec la BAF attendue pour
f = 0,35 — candidat de délétion 22q13 en mosaïque (anneau 22 ?), à vérifier.

Seuil abaissé à 3 fenêtres (3 Mb) le même jour : Alexandre a vu dans IGV une
délétion 22q13 terminale en mosaïque (47,7 Mb-fin, CN 1,62, BAF conforme à
f = 0,38) que le seuil de 5 Mb manquait, comme le caller principal (aucun
appel : CN trop proche de 2 pour le HMM). La dernière fenêtre partielle du
contig est gardée et la fin bornée à sa longueur. Effet : PMS 1 → 3 porteurs,
tous des mosaïques 22q13 terminales (f 0,35-0,56) ; dataset_paper 3 → 4, le
nouveau étant une DUP 7q11.23 de 3 Mb à CN 2,68 (région riche en segdups,
mosaïque ou duplication constitutionnelle sous-estimée — à voir dans IGV).

Limites : autosomes seulement (ploïdie des X/Y absente des stores 0.3) ;
événements < 3 Mb hors de portée ; les délétions complètes (f ≈ 1) n'ont plus
de sites hétérozygotes et ne sont pas confirmées par la BAF (c'est le travail
du caller principal). À faire : intégrer la passe au pipeline Nextflow, et
faire profiter le caller de la fraction (états fractionnaires ou fusion des
fragments d'un segment mosaïque en un appel unique).

### 1terdecies. Drapeau QC du porteur sur les lignes familles (2026-10-01)

Demande d'Alexandre après le rapport dataset_paper : le score d'impact ne
regarde que le contenu génique, et la famille FAM-03 occupait une
dizaine des 57 lignes ≥ +1,00 (PTEN, FGFR2…) avec des enfants à 478 et 94
candidats de novo. Colonne « QC porteur » + filtre (tous / sans alerte / avec
alerte) : wave, charge de sites PASS ≥ 5× la médiane, garde-fou de novo.

Deux défauts trouvés en le validant sur PMS (49 des 58 lignes ≥ +1,00
d'abord signalées, alors que ce sont les vraies délétions SHANK3) :

1. **Le garde-fou de novo à seuil absolu signalait tout le monde** : > 20
   candidats = 664 des 665 enfants dataset_paper et les 77 enfants PMS, et
   leurs tableaux de novo étaient masqués dans le rapport principal. Seuil
   désormais relatif, comme la charge de porteurs : 5× la médiane de la
   cohorte (dataset_paper médiane 75, p90 367, p99 3173 → 63 enfants, dont
   SMP-03 de la famille FAM-03 ; PMS médiane 100 → 1 enfant), le plancher de 20
   restant le seuil des runs famille (< 20 enfants). Le comptage se fait en
   événements (sites fusionnés par la règle du collapse) ; mon hypothèse
   qu'un grand de novo comptait pour des dizaines de sites était fausse,
   l'effet est marginal (SMP-01 : 116 sites, 114 événements).
2. **Le score de wave est gonflé par les grands CNV** : 54 % des porteurs
   d'une DEL chr22 ≥ 1 Mb sont « wavy » contre 13 % des parents ; chr22 seul
   est à 0,4-0,6 chez eux (énergie basse fréquence d'une marche de
   couverture). Exclure les contigs portant un grand appel propre ne marche
   pas (les 319 échantillons PMS en ont un, souvent artefactuel, et le
   signal des vrais wavy s'effondre). Retenu : sur chaque ligne, le wave du
   porteur est la moyenne de ses scores par contig **sans le contig de la
   ligne** (`row_wave_score`) — 26 → 5 porteurs signalés sur leur DEL chr22 ;
   un génome réellement wavy oscille partout et garde son drapeau.

4ᵉ alerte « grands CNV » (2026-10-02) : au moins 3 événements ≥ 1 Mb avec
appel propre PASS chez le porteur (`big_cnv_counts`, par page familiale). Sur
dataset_paper, la famille FAM-03 n'était signalée par aucune des trois alertes :
ses 13 lignes ≥ 1 Mb sont portées par le père SMP-02 (wave 0,112 recalculé
depuis son store 0.3, charge normale, pas d'enfant), 13 délétions non
transmises sur 5 chromosomes — médiane de la cohorte 0, p99 1, max 13 ; le
seuil 3 signale 2 échantillons sur 2673 (SMP-02, SMP-04). Réarrangement
complexe ou problème d'échantillon : à trancher par IGV.

Résultat PMS : lignes ≥ +1,00 avec appel propre PASS, 51 sans alerte / 7 avec
(FAM-09 et FAM-06 portent deux grands CNV, DEL chr22 + DUP chr8 de
translocation ; FAM-14 wave 0,184). Reste à faire : le score de wave
lui-même (stocké au call) devrait exclure les segments CNV appelés — à
recalibrer contre les labels du wave-detector (docs/wave_qc_validation.md).

### 1duodecies. Deux points relevés au rapport dataset_paper du 2026-09-30

- **Fins d'appels au-delà du chromosome (caller, préexistant) — corrigé dans
  le code, données à reprendre.** 24 sites de `cohort/sites.vcf.gz` finissent
  jusqu'à 1 kb après la fin de leur contig (chrX 156,041,002 pour
  156,040,895 ; chr12 133,276,002 pour 133,275,309) : le dernier bin de 1 kb
  d'un segment de profondeur n'était pas borné à la longueur du contig. 1362
  lignes des pages familles en héritent (identique avant et après `ea8a8e3`).
  Correctifs : `447ad86` borne la fin à l'écriture (`vcfout.clamp_end`,
  partagé par `write_vcf` et `run_merge` ; borne interne `length - 1` car
  htslib écrit END = POS + |SVLEN|, cf. (f) : l'END écrit tombe sur la
  dernière base du contig) ; `04d63af` borne les appels propres lus par le
  report (`own_call_index`) ; `_read_cohort_vcf` borne aussi les sites, pour
  les lignes sans appel propre. Reste : re-merge de dataset_paper et PMS
  (`sites.vcf.gz` et `cohort.vcf.gz` publiés débordent encore). Les VCF
  filtrés par échantillon gardent le débordement jusqu'à un rappel : le
  report les borne à la lecture, tout autre lecteur direct les voit encore.
- **Temps du report.** `own_span` + `_cohort_site_pass` font passer les pages
  familles de dataset_paper de 1 h 10 à 2 h 15 (2670 échantillons, 1206
  familles). Piste : résoudre le statut de site de chaque appel propre une
  seule fois par échantillon (et non par famille), ou joindre appels propres
  et sites par un balayage trié par contig plutôt que par bisect par appel.

### 1decies. Les stores d'évidence publiés détruits par publishDir (2026-09-29)

Au relancement de dataset_paper (filter -> report, correctif grandes DEL), les
2670 tâches VEDETTE_FILTER se sont terminées en 10 min avec exit 0 — et
**2071 des 2670 `samples/<S>/<S>.evidence.zarr` publiés sont devenus illisibles**
(`OSError 40: Too many levels of symbolic links`).

Mécanisme : `VEDETTE_FILTER` reçoit le store en entrée et le **ré-émet en
sortie** (pour l'enchaîner vers genotype). Nextflow stage l'entrée par un lien
symbolique vers le chemin publié, puis `publishDir` republie *toutes* les
sorties — dont le store — vers ce même chemin publié. Résultat : le chemin
publié devient un lien vers la copie work, qui est elle-même un lien vers le
chemin publié. Cycle. Un répertoire ne se hardlinke pas, donc `mode: 'link'`
retombe silencieusement sur le lien symbolique et rien n'avertit.

Aucune donnée perdue : les vrais stores vivent dans les répertoires work de
`VEDETTE_CALL`, seul le pointeur avait sauté. Réparation : repointer chaque
lien cassé vers le store réel (choix par nombre de contigs puis mtime), 2071/2071
réparés, 0 irréparable. Vérifications : les 2670 stores s'ouvrent, l'attribut
`sample` correspond au répertoire pour les 2670 (aucun mélange d'échantillons),
et un `vedette filter` rejoué sur trois échantillons réparés redonne le VCF
publié **au bit près** — le store repointé est bien celui qu'avait lu le filtre.

Correctifs : `publishDir ... saveAs: { fn -> fn.endsWith('.evidence.zarr') ? null : fn }`
dans `modules/vedette_filter.nf`, et `tests/test_workflow_modules.py`, qui échoue
si un module republie une entrée qu'il ré-émet (vérifié en remettant l'ancienne
version : le test tombe).

Risque résiduel, non traité : les stores publiés sont des liens vers l'arbre
`workflows/nextflow/work` (2,6 To au 2026-09-29). Un `nextflow clean` les détruirait, et il faudrait
rappeler depuis les CRAM (~6 jours sur 2670). Deux options le jour où ça compte :
publier les stores en `copy` (coût disque), ou ne pas les publier du tout et
assumer le work dir comme emplacement canonique.

### 1nonies. Les grandes délétions tombaient sous le seuil SQ (2026-09-28)

Alexandre, sur FAM-01-001 : IGV montre une délétion terminale nette
chr22:41,9-50,8 Mb chez le proposant, Vedette n'en rapportait que des morceaux.

Les appels propres couvraient pourtant toute la région — mais avec des SQ de
3,73 et 4,31 sur les deux plus gros segments (6,43 Mb et 229 kb), donc
`LowSQ`, donc absents du merge. Mesure sur 25 génomes PMS : **50 % des DEL de
1 à 10 Mb ont SQ < 5**.

Cause, en deux temps :
1. Le modèle n'a presque pas appris à ces tailles : 131 exemples > 1 Mb sur
   132 555 DEL étiquetées (0,10 %), et **zéro** exemple DUP vrai au-dessus de
   100 kb. Pire, ces 131 ont des features dégénérées (dhffc = 0, fc_chrom = 0,
   CN non appelé) : ce sont des zones de perte totale de couverture, pas des
   délétions hétérozygotes.
2. Quand le HMM fractionne une grande délétion en segments adjacents, **les
   flancs de chaque segment tombent dans le segment voisin**, donc
   `dhffc` ≈ 1 : la feature la plus discriminante du modèle devient aveugle.
   Sur le cas réel, dhffc = 1,03 alors que `fc_chrom` = 0,528 dit clairement
   « une copie ».

Correctif : une DEL d'au moins `large_del_rescue_size` dont le `fc_chrom`
confirme la perte (< 0,75) et peu masquée garde son PASS quel que soit le SQ,
avec `INFO/SQ_OVERRIDE=LARGE_DEL_DEPTH` pour la traçabilité. Seuil calé sur
les données étiquetées : à 200 kb le rattrapage prend 5 appels, tous vrais ;
à 100 kb il en prend 23 dont 10 faux ; à 1 Mb il rate le premier segment de
229 kb des délétions terminales PMS. Sur SMP-07, 2 rattrapages sur le génome
entier et l'étendue PASS passe de fragments épars à 41 894 001-50 799 001.

À faire : c'est un garde-fou, pas une correction de fond. Il faut (a) des
exemples de grands événements réalistes à l'entraînement v6 — les délétions
terminales de PMS en fournissent, à condition de les valider indépendamment ;
(b) corriger la contamination des flancs de `dhffc` (exclure des flancs les
bins couverts par un autre appel du même échantillon) ; (c) rien n'est fait
pour les grandes DUP, faute du moindre exemple vrai à l'entraînement.

### 2bis. Score d'impact plus discriminant (demande Alexandre, 2026-09-28) — v1 FAIT (2026-09-30)

**v1 en place** (`_GeneIndex.impact`, `_gene_points`, `_GeneTrek`) : points de
contenu génique d'après la grille ACMG/ClinGen des CNV (Riggs et al. 2020,
sections 1-3), séparés DEL / DUP, GeneTrek tenant lieu de carte de dosage
ClinGen. Gène « type haploinsuffisant » = NDD établi (`hc_ndd_v7` ou
`sysndd_definitive`) à transmission dominante (l'une des trois colonnes
`inheritance*`). DEL : complète 1,00 / exonique 0,90 ; NDD établi mode inconnu
0,45 ; récessif seul 0,15 ; NDD candidat (candidate_ndd, SFARI 1, autisme,
PanelApp ID/épilepsie/autisme) 0,30, 0,45 si contraint ; contraint seul
(`is_shet_constrained`, s_het ≥ 0,1 ou LOEUF < 0,35) 0,15 ; section 3 :
≥ 25 / ≥ 35 gènes codants +0,45 / +0,90. DUP : DUP intragénique d'un gène NDD
dominant 0,45 (2I, perte de fonction possible), DUP complète 0,15 au plus
(triplosensibilité absente de GeneTrek) ; ≥ 35 / ≥ 50 gènes +0,45 / +0,90.
Aucun gène codant −0,60 (1B). Libellé : points, niveau, règle, gène moteur
avec ses preuves (NDD établi/candidat, AD/AR, s_het), autres gènes NDD, nombre
de gènes codants. Tri numérique dans les pages familles (`data-impact`).

Mesuré sur les 576 429 lignes PASS des pages familles PMS : l'ancien score
avait 5 valeurs et **3 664 lignes au niveau maximal** ; le nouveau en a 10,
et **49 lignes à +1,90 — exactement les délétions 22q13 terminales emportant
SHANK3 avec ≥ 35 gènes codants**. Suivent +1,45 / +1,00 (délétions SHANK3 plus
courtes), +1,05 (duplications 8q24.3 de 75-103 gènes des translocations
déséquilibrées, FAM-09 et FAM-06), puis 499 atteintes exoniques de gènes NDD
dominants (+0,90, beaucoup sans appel propre).

Reste pour une v2 : triplosensibilité (ClinGen dosage map, à ajouter au
fichier), 5'/3' du transcrit MANE pour les DEL partielles (2C vs 2D-2E), et
sections 4-5 (de novo, ségrégation, phénotype de la cohorte) — à caler sur les
de novo déjà validés par IGV.

**Formulation initiale de la demande :**


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

### 2ter. Vérification au niveau du génotype — FAIT, mais pas à la duphold (2026-09-28)

Testé avant d'implémenter, sur les 77 trios PMS avec le taux de de novo
apparent comme proxy d'erreur de génotypage. **Deux approches duphold ont été
mesurées et écartées** :

- **DHFFC (fold-change contre les flancs)**, le cœur de duphold : ne
  discrimine rien ici, il est même inversé pour les DEL (1,4 % d'erreur quand
  il est « cohérent » contre 0,4 % sinon). La raison : **98 % des porteurs DUP
  sont à des sites de moins de 3 kb**, où la profondeur n'a pas de résolution
  (bins de 1 kb) — le fold-change y est du bruit.
- **Exiger des clips des deux côtés** du point de cassure : inversé aussi
  (DUP 10,5 % d'erreur à deux côtés contre 5,8 % à un seul).

Ce qui marche, mesuré :

| critère (parmi porteurs GQ>=20) | DEL gardés / erreur | DUP gardés / erreur |
|---|---|---|
| aucun | 100 % / 0,8 % | 100 % / 1,9 % |
| **CN cohérent avec le type** | 88 % / **0,3 %** | 24 % / **0,3 %** |
| RD cohérent | 92 % / 0,6 % | 26 % / 0,6 % |
| CN cohérent ou PE>=2 | 92 % / 0,4 % | 68 % / 1,3 % |

Le **CN re-génotypé** bat le RD et ne coûte rien (déjà chargé par le rapport).
Le canal jonction est confirmé peu spécifique : « SR seul » est la pire
catégorie (petites DEL 5,9 % d'erreur contre 0,8 % avec PE>=2 ; petites DUP
6,9 % contre 3,0 %).

Fait : colonne **support** (« fort » = le CN d'un porteur confirme le type)
sur les pages events, avec sélecteur. **Pas** un filtre par défaut : sur les
42 sites relus sous IGV, l'exiger fait tomber 2 des 3 vrais (précision 23 %
-> 20 %). Les deux mesures se contredisent sur ce point — le proxy trio dit
que c'est un bon indicateur d'erreur de génotypage, l'échantillon IGV (3 vrais,
donc très fragile) dit qu'il écarte du vrai. À retrancher quand il y aura plus
de validations.

Reste à faire : rendre `_junction_support` spécifique (apparier les clips à la
jonction — le store ne garde pas le partenaire SA, donc cela demande de
l'ajouter à l'extraction), et un DHBFC (bins de GC comparable) si un jour la
profondeur redevient le canal limitant sur de gros événements.

### 2quater. (archive) formulation initiale de la demande duphold

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
