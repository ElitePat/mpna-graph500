#!/bin/bash

# ============================
# Configuration
# ============================
PROGRAM=./graph500_custom_bfs
SCALE=20
EDGE_FACTOR=16

# Liste de cœurs à tester
CORES_LIST=(1 2 4 8 16 32 64 128)

# Dossier pour stocker les logs
LOG_DIR=logs4

mkdir -p $LOG_DIR

# ============================
# Boucle sur chaque nombre de cœurs
# ============================
for NP in "${CORES_LIST[@]}"
do
    LOG_FILE="$LOG_DIR/run_${NP}cores.log"
    echo "=== Lancement sur $NP cœurs ==="
    
    # Lancer le programme avec mpirun et rediriger la sortie vers un log
    mpirun -np $NP $PROGRAM $SCALE $EDGE_FACTOR > $LOG_FILE 2>&1

    # Vérifier si le run s'est bien terminé
    if [ $? -eq 0 ]; then
        echo "Run terminé avec succès : $LOG_FILE"
    else
        echo "Erreur lors du run sur $NP cœurs"
    fi
done

echo "Tous les runs sont terminés. Les logs sont dans $LOG_DIR"
