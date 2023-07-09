#!/bin/bash
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 0b9606e (dash app container)
=======
>>>>>>> 0b9606e (dash app container)
source activate mod-site


python ./dash_main_n.py
# gunicorn -w 6 --threads=12 --worker-class=gthread -b 0.0.0.0:5000 --timeout 120 --max-requests 500 --max-requests-jitter 100 --graceful-timeout 120 dash_main_n:server --access-logfile /app/logs/access.log
<<<<<<< HEAD
<<<<<<< HEAD
=======
conda activate mod-site && gunicorn -w 6 --threads=12 --worker-class=gthread -b 0.0.0.0:5000 --timeout 120 --max-requests 500 --max-requests-jitter 100 --graceful-timeout 120 app:server --access-logfile /app/logs/access.log
>>>>>>> 0cfd125 (added docker)
=======
>>>>>>> 0b9606e (dash app container)
=======
>>>>>>> 0b9606e (dash app container)
