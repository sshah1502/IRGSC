from astroquery.ukidss import Ukidss
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii
from astropy.table import Table
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import json
import requests
try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve
try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib 
import mastcasjobs

username='sshah1502'
password = 'Astroboy1502'

class ukidss(object):
    def __init__(self, ra, dec, search_radius):
        """
        Function to fetch UKIDSS observed NIR data to validate the computed NIR magnitudes
        Inputs required as ra, dec in degrees and search_radius in arcmin.
        The frame of reference is 'icrs' which can be changed if required.
        UKIDSS has several subsurveys like:
        Large Area Survey - LAS
        Galactic Plane Survey - GPS
        Galactic Cluster Survey - GCS
        Extragalactic Deep Survey - DXS
        Ultra Deep Survey - UDS
        The argument programme_id requires the name of one these subsurveys and .
        """
        sub_surveys = ['LAS', 'GPS', 'GCS', 'DXS', 'UDS']
        for i in range(len(sub_surveys)):
            try:
                table = Ukidss.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg),\
                                                       frame='galactic'),\
                                                    radius=search_radius*u.arcmin, programme_id = i)
            except NoResultsWarning:
                print("No results")

class panstarrs():
  def __init__(self, ra, dec, search_radius):
    #if self.get_optical_data is True:
          self.ra,self.dec, self.search_radius = ra, dec, search_radius

          self.jobs = mastcasjobs.MastCasJobs(username=username, password=password, userid=None,\
                request_type="GET", context="PanSTARRS_DR2", base_url="https://mastweb.stsci.edu/ps1casjobs/services/jobs.asmx",\
                    wsid_url=None, fast_url=None)
          if 'mydb' in self.jobs.list_tables():
           self.jobs.drop_table('mydb')

          self.query = f"""
                            select o.objid, o.raMean, o.decMean, o.raMeanErr, o.decMeanErr,
                            m.gPSFMag, m.gPSFMagErr, m.gKronMag, m.gKronMagErr, m.rPSFMag,
                            m.rPSFMagErr, m.rKronMag, m.rKronMagErr, m.iPSFMag, m.iPSFMagErr,
                            m.iKronMag, m.iKronMagErr, m.zPSFMag, m.zPSFMagErr, m.zKronMag,
                            m.zKronMagErr, m.yPSFMag, m.yPSFMagErr, m.yKronMag, m.yKronMagErr,
                            o.objInfoFlag, o.qualityFlag, o.nDetections, o.nStackDetections,
                            m.ginfoFlag, m.ginfoFlag2, m.ginfoFlag3, m.rinfoFlag, m.rinfoFlag2,
                            m.rinfoFlag3, m.iinfoFlag, m.iinfoFlag2, m.iinfoFlag3, m.zinfoFlag,
                            m.zinfoFlag2, m.zinfoFlag3, m.yinfoFlag,
                            m.yinfoFlag2, m.yinfoFlag3

                            into mydb

                            from fGetNearbyObjEq({self.ra}, {self.dec}, {self.search_radius}) nb
                            join ObjectThin o on o.objid=nb.objid
                            join StackObjectThin m on o.objid=m.objid
                            join StackObjectView s on s.objid=o.objid and
                            o.uniquePspsOBid = s.uniquePspsOBid
                            join AstrometryCorrection a on a.objid=o.objid
                            WHERE s.bestDetection >0 and s.primaryDetection > 0
                         
          """
          self.id = self.jobs.submit(self.query, task_name="t0")
          self.status = self.jobs.status(self.id)
          print(self.status, self.id)


  def get_output(self):
          if self.status[1] == 'finished':
              print(self.status[1])
              print(self.jobs.list_tables(context="MYDB"))
              self.jobs.request_and_get_output('t0','CSV'
                                               ,'PS1' + '_' + str(self.ra) + '_' + str(self.dec) + '.csv')
              df = pd.read_csv('PS1' + '_' + str(self.ra) + '_' + str(self.dec) + '.csv')
              self.df = df

  def mastQuery(self,request, json_return=False):
        """
        Perform a MAST query.

        Parameters
        ----------
        request (dictionary): The MAST request json object

        Returns the text response or (if json_return=True) the json response
        """

        url = "https://mast.stsci.edu/api/v0/invoke"

        # Encoding the request as a json string
        requestString = json.dumps(request)

        # make the query
        r = requests.post(url, data=dict(request=requestString))

        # raise exception on error
        r.raise_for_status()

        if json_return:
            return r.json()
        else:
            return r.text


  def resolve(self,name):
      """Get the RA and Dec for an object using the MAST name resolver

      Parameters
      ----------
      name (str): Name of object

      Returns RA, Dec tuple with position"""
      resolverRequest = {'service':'Mast.Name.Lookup',
                      'params':{'input':name,
                                'format':'json'
                                },
                      }
      resolvedObject = self.mastQuery(resolverRequest, json_return=True)
    # The resolver returns a variety of information about the resolved object, 
    # however for our purposes all we need are the RA and Dec
      try:
          objRa = resolvedObject['resolvedCoordinate'][0]['ra']
          objDec = resolvedObject['resolvedCoordinate'][0]['decl']
      except IndexError as e:
          raise ValueError("Unknown object '{}'".format(name))
      return (objRa, objDec)
      
     

ps = panstarrs(0,0,30)
ps.get_output()