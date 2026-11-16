# Wave QC — validation du score (2026-09-14)

`svengers.qc.sample_wave_score` mesure la fraction d'énergie spectrale basse
fréquence (périodes >= 200 kb, fenêtres de 10 Mb, détrend linéaire) du ratio de
couverture **corrigé GC** — le signal qui alimente réellement le HMM. Vérité
externe : classifications du pipeline GHFC `Wave_hg38_remapping_slurm_nextflow`
(seuil = 95e percentile de calibration sur 4017 bedgraphs bruts).

## Protocole

592 génomes de la cohorte dataset_paper avec un evidence store publié :
les 172 étiquetés WAVES disponibles + 420 NORMAL tirés au hasard (seed 1).
Score calculé sur `depth_corr / mean_depth` par contig (bins 1 kb, masques
exclude/N appliqués), moyenne des médianes par fenêtre.

## Résultats

| métrique | valeur |
|---|---|
| AUC vs labels GHFC | **0,946** |
| Spearman vs score composite GHFC | 0,651 |
| médiane NORMAL / WAVES | 0,112 / 0,174 |
| p95 NORMAL | 0,145 |

Seuils candidats (sur cet échantillon, enrichi à 29 % de waves) :

| seuil | sensibilité | spécificité |
|---|---|---|
| 0,128 (Youden) | 0,93 | 0,87 |
| **0,145 (p95 normaux — retenu, `qc.WAVE_WARN`)** | **0,80** | **0,95** |
| 0,160 | 0,67 | 0,97 |

La dynamique est compressée par rapport aux scores du détecteur GHFC (calculés
sur la couverture brute) : la correction GC absorbe une partie de l'artefact.
Le lien opérationnel avec le runtime est documenté dans la même analyse :
les échantillons WAVES ont une durée de `svengers call` médiane 5× supérieure
(8,5 h vs 1,7 h, odds ratio 8,8 pour le décile le plus lent), d'où le budget
`breakpoint.resolve_max_alignments` introduit en même temps que ce score.

Reproduire : `scratch validate_wave.py` (session), stores dataset_paper +
`final_results/merged_summary.csv` du pipeline Wave.

## Mise à jour 2026-09-25 : seuil abaissé à 0,128 (Youden)

Le garde-fou de novo (report.py) a laissé passer deux échantillons wavy à
0,139 et 0,140 (SMP-06/X), sous le seuil p95 de 0,145, avec le même profil
qu'un cas confirmé (SMP-05, 0,178) — génotypage globalement peu fiable,
des milliers de faux « de novo » chez leur enfant malgré tous les filtres de
site. Le seuil `qc.WAVE_WARN` passe au point de Youden déjà documenté plus
haut (0,128 ; sensibilité 0,93, spécificité 0,87), au prix d'un peu plus de
faux positifs sur les échantillons réellement normaux.
