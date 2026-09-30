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

### 1duodecies. Deux points relevés au rapport dataset_paper du 2026-09-30

- **Fins d'appels au-delà du chromosome (caller, préexistant).** 24 sites de
  `cohort/sites.vcf.gz` finissent jusqu'à 1 kb après la fin de leur contig
  (chrX 156,041,002 pour 156,040,895 ; chr12 133,276,002 pour 133,275,309) :
  le dernier bin de 1 kb d'un segment de profondeur n'est pas borné à la
  longueur du contig. 1362 lignes des pages familles en héritent (identique
  avant et après `ea8a8e3`). À corriger à l'écriture du VCF (vcfout : borner
  `stop` à la longueur du contig), puis re-merge ; le report peut borner en
  attendant.
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
