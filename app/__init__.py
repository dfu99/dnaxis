from flask import Flask
from flask_session import Session
from flask_mail import Mail

from app import config


def create_app():
    _app = Flask(__name__)

    # Mail config
    _app.config['MAIL_SERVER'] = 'smtp.gmail.com'
    _app.config['MAIL_PORT'] = 465
    _app.config['MAIL_USERNAME'] = 'dnaxis.webmaster@gmail.com'
    _app.config['MAIL_PASSWORD'] = config.MAIL_PASSWORD
    _app.config['MAIL_USE_TLS'] = False
    _app.config['MAIL_USE_SSL'] = True
    Mail(_app)

    # Session config
    _app.config['SESSION_TYPE'] = 'redis'
    Session(_app)

    # Upload directory
    _app.config['UPLOAD_FOLDER'] = 'jobs'

    # Secret key for flash messages
    _app.secret_key = b'kxxtMg!5BG3&Z9Rp'

    # Register blueprint
    from app.routes import bp
    _app.register_blueprint(bp)

    return _app


# Module-level app for backward compatibility
app = create_app()
