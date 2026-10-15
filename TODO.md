# TODO

- [ ] Revoir la gestion des échecs de BRAKER2, BRAKER3 et GALBA, notamment ceux liés à Augustus. Pour l'instant, conserver le comportement actuel : après cinq tentatives infructueuses, produire un GFF vide, écarter cette prédiction et laisser le pipeline continuer avec les autres méthodes. Prévoir ultérieurement un bilan explicite des méthodes en échec pour distinguer une exécution complète d'une exécution avec des prédictions manquantes.

- [ ] Valider en exécution les nouvelles étapes de sélection automatique, annotation non codante et export final. Pour cette étape, seule une revue statique a été effectuée à la demande de l'utilisateur. Vérifier notamment les fichiers Rfam/OMA et les rapports de validation du dépôt cible avant tout dépôt.

- [ ] Adapter ultérieurement le cache aux changements d'options des programmes et de banques. Ne pas modifier son comportement pour cette étape, conformément à la demande de l'utilisateur.
