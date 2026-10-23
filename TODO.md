# TODO

- [ ] Revoir la gestion des échecs de BRAKER2, BRAKER3 et GALBA, notamment ceux liés à GeneMark ou Augustus. Pour l'instant, conserver le comportement actuel : après cinq tentatives infructueuses, produire un GFF vide, écarter cette prédiction et laisser le pipeline continuer avec les autres méthodes. Prévoir ultérieurement un bilan explicite des méthodes en échec pour distinguer une exécution complète d'une exécution avec des prédictions manquantes.

- [x] Valider la sélection automatique, l’annotation non codante et l’export final sur le plus petit profil `test` (Hanseniaspora uvarum, trois chromosomes). Exécution complète et reprise réussies le 23 octobre 2026 ; voir `test/validation/minimal_validation_20261023.txt`.

- [ ] Étendre la validation aux trois profils minimaux par genre, aux profils full et aux comparaisons entre plusieurs génomes. Effectuer la validation avec les outils des banques et les métadonnées enregistrées avant tout dépôt.

- [ ] Adapter ultérieurement le cache aux changements d'options des programmes et de banques. Ne pas modifier son comportement pour cette étape, conformément à la demande de l'utilisateur.
