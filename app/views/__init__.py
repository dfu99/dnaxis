from app import app
from flask_session import Session
from flask_mail import Mail
from app import config
"""Mail config"""
app.config['MAIL_SERVER'] = 'smtp.gmail.com'
app.config['MAIL_PORT'] = 465
app.config['MAIL_USERNAME'] = 'dnaxis.webmaster@gmail.com'
app.config['MAIL_PASSWORD'] = config.MAIL_PASSWORD
app.config['MAIL_USE_TLS'] = False
app.config['MAIL_USE_SSL'] = True

mail = Mail(app)

SESSION_TYPE = 'redis'
app.config.from_object(__name__)
Session(app)

# set where to save uploads
UPLOAD_DIR = 'jobs'
app.config['UPLOAD_FOLDER'] = UPLOAD_DIR

# default scaffolds and their lengths
scafLen = {'m13mp18': 7249, 'phix174': 5386, 'p8064': 8064, 'p7308': 7308, 'debug':50000}

# for flash requests
app.secret_key = b'kxxtMg!5BG3&Z9Rp'

from .connections import *
from .info import *
from .pathway import *
from .stl import *
from .draw import *
from .results import *
from .submit import *
from .verify import *