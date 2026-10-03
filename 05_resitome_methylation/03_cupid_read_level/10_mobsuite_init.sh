#!/bin/bash
set -euo pipefail
. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate cupid_mobsuite
mob_init -v
