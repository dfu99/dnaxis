# DNAxiS

DNAxiS (caddna.cs.duke.edu) is a browser-based app for DNA origami design of shapes with axial symmetry.

The server is written using Python 3.7 and Flask.

## To run the webserver locally:
### Prerequisites
(For Windows users: We recommend using Git Bash or another Linux-based terminal of your choice)

Install Redis (https://redis.io) and run redis-server.exe

### Setup
Create a virtual environment in the local directory:
  python -m venv <your_venv_name>
 
Then activate it.

Windows: 

    source <your_venv_name>/Scripts/activate
Mac/Linux:

    source <your_venv_name>/bin/Activate

Install the packages as detailed in the requirements.txt

    pip install -r requirements.txt

Then run the server using

    python DNAxiS-websrv.py

You can access it locally by opening a web browser and going to the URL 127.0.0.1:5000


