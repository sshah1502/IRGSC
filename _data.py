import pandas as pd
import query
import mastcasjobs
class casjobs(object):
    def __init__(self,name = None, radius = 30):
        if name is not None:
          self.name        = name
          self.ra,self.dec = self.resolve(self.name)

          self.jobs = mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
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

                            from fGetNearbyObjEq({self.ra}, {self.dec}, {radius}) nb
                            join ObjectThin o on o.objid=nb.objid
                            join StackObjectThin m on o.objid=m.objid
                            join StackObjectView s on s.objid=o.objid and
                            o.uniquePspsOBid = s.uniquePspsOBid
                            join AstrometryCorrection a on a.objid=o.objid
                            WHERE 
                            s.bestDetection > 0 and 
                            s.primaryDetection > 0 and
                            mag.gPSFMag!=-999 and 
                            mag.rPSFMag!=-999 and 
                            mag.iPSFMag!=-999 and 
                            mag.zPSFMag!=-999 and 
                            mag.yPSFMag!=-999 and
                            mag.gKronMag!=-999 and 
                            mag.rKronMag!=-999 and 
                            mag.iKronMag!=-999 and 
                            mag.zKronMag!=-999 and 
                            mag.yKronMag!=-999 and
                            mag.gPSFMagErr!=-999 and
                            mag.rPSFMagErr!=-999 and
                            mag.iPSFMagErr!=-999 and
                            mag.zPSFMagErr!=-999 and
                            mag.yPSFMagErr!=-999 and

                            mag.gPSFMag - mag.gKronMag <0.05 and 
                            mag.rPSFMag - mag.rKronMag <0.05 and
                            mag.iPSFMag - mag.iKronMag <0.05 and
                            mag.zPSFMag - mag.zKronMag <0.05 and
                            mag.yPSFMag - mag.yKronMag <0.05
                         
          """
          self.id = self.jobs.submit(self.query, task_name="python cone search")
          self.status = self.jobs.status(self.id)
    def get_output(self):
      self.status  = self.jobs.status(self.id)
      if self.status[1]== 'finished':

        self.jobs.request_and_get_output('mydb','CSV'
        ,'/content/INSIST/data/output.csv')

        df = pd.read_csv('/content/INSIST/data/output.csv')
        df = df.rename(columns = {'raStack' : 'ra', 'decStack' : 'dec'})

        for i in ['g','r','i','z','y']:
          continue

          df[f'{i}Flux']     = 3631*pow(10,-df[f'{i}PSFMag']/2.5)*1000             # mJy
          df[f'{i}Flux_err'] = (df[f'{i}PSFMagErr']*df[f'{i}Flux'])/1.082        # mJy

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
      
    def perform_mag_cut(self, threshold_error):

      for m_err in ['g', 'r',	'i',	'z', 'y']:
        self.df = self.df[self.df[f'{m_err}PSFMagErr']]