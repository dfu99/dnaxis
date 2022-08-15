from flask import Flask

app = Flask(__name__)

from app import config
from app import views
from app import forms