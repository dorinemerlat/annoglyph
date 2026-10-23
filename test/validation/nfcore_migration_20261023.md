# Migration complémentaire vers nf-core — 23 octobre 2026

État : tests ciblés terminés avec succès sous Nextflow 25.10.2,
Singularity et Slurm, avec `conf/ifb.config` et
`test/validation/slurm_runtime.config`.

## Changements

Barrnap, Infernal/cmscan, tRNAscan-SE, BUSCO protéines, STAR et la récupération
SRA utilisent maintenant les modules nf-core. Les nouveaux modules sont figés
sur `ffd1a93ab6f47daae9b0b834ecf024c7af6e223e`. BUSCO réutilise le module déjà
installé dans le dépôt. Révisions et adaptations sont suivies dans `modules.json`.
BRAKER2, BRAKER3 et GALBA conservent leur implémentation locale.

Les images des logiciels existantes sont conservées par la configuration pour
éviter de mêler migration des modules et changement de version des outils.
Les options biologiques restent les mêmes. Les adaptations concernent les
sorties compressées d’Infernal, les structures et le journal optionnel de
tRNAscan-SE, les versions des outils et le chemin de banque BUSCO.

## Vérifications

- Test d’ensemble des modules migrés : 13 tâches terminées ou reprises du cache,
  exécution finale `golden_angela`, session
  `42d0881c-c66b-477b-934e-4511745bc417`, sortie 0.
- Génome : extrait Hanseniaspora uvarum de 2 601 250 bases, trois contigs.
- ncRNA : 30 tRNA et 6 autres ARN non codants, soit 36 annotations biologiques.
  Les 70 lignes de caractéristiques GFF, dont 34 exons, sont strictement
  identiques aux résultats du module local précédent.
- tRNAscan-SE : GFF, FASTA compressé, statistiques et structures produits.
  Avec `--progress`, la version 2.0.12 n’écrit pas le journal demandé avec `-l` :
  sa sortie est donc optionnelle. Le FASTA compressé est correctement déclaré.
- BUSCO protéines : résumé JSON produit et ensemble du champ `results`
  identique au test précédent : 202 BUSCO complets, 6 fragmentés, 2 111 absents,
  sur 2 319 marqueurs. Banque locale saccharomycetes_odb12, mode hors ligne.
- Le dossier de téléchargement BUSCO est conservé avec sa hiérarchie de
  lignées : renommer directement une lignée en `input_lineages` faisait chercher
  `input_lineages_odb12` à BUSCO 6. Les échecs de lot en mode protéines sont
  maintenant des erreurs explicites, même si BUSCO renvoie initialement 0.
- STAR : indexation, alignements simples et appariés, BAM triés et journaux
  vérifiés sur 5 000 lectures/spots par jeu. Index de suffixes fixé à 9 pour
  ce petit génome, comme auparavant. La collecte des versions ne demande pas
  GNU Awk, absent de l’image existante et inutile avec cet index explicite.
- Téléchargement et conversion SRA : petit run simple DRR000774 récupéré ;
  conversion simple et conversion appariée vérifiées. Le cas apparié utilise
  l’archive réduite SRR11140744 des données de test nf-core. Les métadonnées
  du runlist déterminent le type de lectures dans le pipeline actif.
- Sous-workflow RNA-seq complet : dix tâches réussies avec deux identifiants
  du même assemblage. `sample_a` reçoit un jeu simple et un jeu apparié,
  `sample_b` un jeu simple. Deux index sont créés, puis trois alignements,
  et les BAM de sample_a sont fusionnés. Les sorties correspondent bien
  aux deux identifiants ; 2 722 et 1 190 alignements sont comptés respectivement.
- Huit tests Python existants passent : sous-échantillonnage (5) et parsing des
  métadonnées SRA (3).
- Chargement du workflow principal avec le profil `singularity,test` et les
  deux configurations IFB : succès en mode preview.
- Les neuf patches nf-core concernés passent `git apply --reverse --check`.
  `git diff --check` est propre.

## Traces et limites

Les journaux et traces ciblés sont dans `.apptainer_tmp/nfcore_migration/` :
`final_run.log`, `trace.tsv`, `busco.log`, `trace_busco.tsv`,
`rna_workflow.log`, `trace_rna_workflow.tsv`, `trna_sra.log` et
`trace_trna_sra.tsv`. Les fichiers de validation sont isolés des résultats
publiés du précédent test minimal. La reprise finale utilise un UUID explicite
et un historique de validation séparé des autres exécutions du projet.

Il s’agit de validations des modules remplacés et du sous-workflow RNA-seq,
avec comparaison des résultats ncRNA et BUSCO. Le pipeline complet, les
prédicteurs de gènes et les jeux full n’ont pas été relancés pour cette migration.
Les sorties des anciennes implémentations restent intactes ; les outils migrés
ont un espace de cache distinct `nfcore_tools_v2`.
