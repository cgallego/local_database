# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:10:52 2014

@ author (C) Cristina Gallego, University of Toronto
"""

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


# configure Session class with desired options
Session = sessionmaker()
myengine = create_engine('sqlite:///C:\\Users\\windows\\Documents\\repoCode-local\\local_database\\localData.db', echo=False)

# later, we create the engine
Base = declarative_base(myengine)