'''
Created on Oct 31, 2014
Updated Mar 2022

@author: omarcillo, fkdd, jwebster
'''

from infrapy.database.taskbase import Base

import sys

from pisces.util import load_config
from obspy.core import UTCDateTime
from datetime import datetime

import numpy as np
from IPython import embed

from sqlalchemy import func

from infrapy.utils.short_time import short_time
from infrapy.association import hjl
from infrapy.utils import data_io
from infrapy import schema

import numpy as np

import pathos.multiprocessing as mp
from multiprocessing import cpu_count

import warnings
from sqlalchemy import exc as sa_exc

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=sa_exc.SAWarning)

class AssocInfraPy_LANL(Base):
    '''
    classdocs
    '''

    algorithm='Blom and Euler'

    def __init__(self, conf_file=[]):
        '''
        Constructor
        '''
        super(AssocInfraPy_LANL,self).__init__(conf_file,'AssocLocParams')
        self.assocversion=0
        print('Assoc version:',self.assocversion)
        self.year=int(self.general_PARAM['year'])
        self.dayofyearini=int(self.general_PARAM['dayofyearini'])
        self.dayofyearend=int(self.general_PARAM['dayofyearend'])
        self.jdayini = int(self.year*1000 + self.dayofyearini)
        self.jdayend = int(self.year*1000 + self.dayofyearend)
        self.cpu=self.general_PARAM['cpucnt']
        self.pl = mp.ProcessingPool(cpu_count() - 1)


        self.net=self.task_PARAM['network']
        self.pfdid=self.task_PARAM['pfdetectid']
        self.pfkid=self.task_PARAM['pfkid']
        self.beamwidth=float(self.task_PARAM['beamwidth'])
        self.rangemax=float(self.task_PARAM['rangemax'])
        self.distmax=float(self.task_PARAM['distmax'])
        self.clusterthresh=float(self.task_PARAM['clusterthresh'])
        self.trimthresh=(self.task_PARAM['trimthresh'])
        self.trimthreshscalar=(self.task_PARAM['trimthreshscalar'])
        self.mindetpop=(self.task_PARAM['mindetpop'])
        self.minarraypop=(self.task_PARAM['minarraypop'])
        self.duration=float(self.task_PARAM['duration'])


        try:
            self.resultstable=self.task_PARAM['resultstable']
        except :
            print('No specific tables')
            self.resultstable=None

        listK=list(self.task_PARAM.keys())
        self.fdtables_names=[]
        for li in listK:
            if bool('fdtable_' in li):
                self.fdtables_names.append(self.task_PARAM[li])
        if len( self.fdtables_names)==0:
            print('NO tables with fd results were included, define one fd results tables, or include Fd_results (this is the table where results are written by default, but needs to be specified ) ')
            sys.exit()
        self.num_tables=len(self.fdtables_names)
            #self.fdtables.append('Fd_results')
        #embed()

    def database_connecting(self):
        print('connecting')

        session,tables = load_config(self.db_PARAM)
        self.session=session
        self.Site=tables['site']
        self.Wfdisc=tables['wfdisc']
        self.Affiliation=tables['affiliation']
        import pisces.schema.css3 as kba

        #class FK_results(schema.fk_results):
        #    __tablename__ = 'FK_results'

        class FK_params(schema.fk_params):
            __tablename__ = 'FK_params'

        class FD_params(schema.fd_params):
            __tablename__ = 'FD_params'

        self.dict_namefk={}

        self.fdtables=[]
        for  fdi in self.fdtables_names:
            self.fdtables.append(type(str(fdi),(schema.fd_results,),{'__tablename__':str(fdi)}))


        class ASSOC_params(schema.ASSOC_params):
            __tablename__ = 'ASSOC_params'

        class ASSOC_results(schema.ASSOC_results):
            __tablename__= self.resultstable


        self.FK_par=FK_params
        #self.FK_results=FK_results
        self.FD_par=FD_params
        #self.Fd_results=Fd_results

        self.ASSOC_par=ASSOC_params
        self.ASSOC_results=ASSOC_results

        self.ASSOC_par.__table__.create(self.session.bind,checkfirst=True)
        self.ASSOC_results.__table__.create(self.session.bind,checkfirst=True)


        try:
            self.Passoc_Q=self.session.query(self.ASSOC_par). \
                filter(self.ASSOC_par.beamwidth==self.beamwidth).\
                filter(self.ASSOC_par.rangemax==self.rangemax). \
                filter(self.ASSOC_par.clusterthresh==self.clusterthresh).\
                filter(self.ASSOC_par.trimthresh==self.trimthresh).\
                filter(self.ASSOC_par.trimthreshscalar==self.trimthreshscalar).\
                filter(self.ASSOC_par.mindetpop==self.mindetpop).\
                filter(self.ASSOC_par.minarraypop==self.minarraypop).\
                filter(self.ASSOC_par.duration==self.duration).\
                all()


            if len(self.Passoc_Q)>1:
                print('issue with the database too many parameters entries, there should be just one')
                embed()
            if len(self.Passoc_Q)==1:
                self.Passoc_Q=self.Passoc_Q[0]
        except Exception as x1:
            print("issue with the table or first assoc entered")
            print(x1)
            embed()
            self.Passoc_Q=[]
            print(self.Passoc_Q)
        if bool(self.Passoc_Q)==False:
            print('New process parameters, write process to INFRA_ASSOC_PARAM table')
            new_row=self.session.query(self.ASSOC_par).count()
            try:
                res=self.ASSOC_par(  beamwidth=self.beamwidth,\
                                     rangemax=self.rangemax,\
                                     clusterthresh=self.clusterthresh,\
                                     trimthresh=self.trimthresh,\
                                     trimthreshscalar=self.trimthreshscalar,\
                                     mindetpop=self.mindetpop,\
                                     minarraypop=self.minarraypop,\
                                     algorithm=self.algorithm,\
                                     duration=self.duration,\
                                     passocid=new_row)
            except Exception as x1:
                print('problem writing to the assoc param file')
                print("Unexpected error:", x1)
                embed()

            self.session.add(res)
            self.session.commit()
            self.Passoc_Q=self.session.query(self.ASSOC_par). \
                filter(self.ASSOC_par.beamwidth==self.beamwidth).\
                filter(self.ASSOC_par.rangemax==self.rangemax). \
                filter(self.ASSOC_par.clusterthresh==self.clusterthresh).\
                filter(self.ASSOC_par.trimthresh==self.trimthresh).\
                filter(self.ASSOC_par.trimthreshscalar==self.trimthreshscalar).\
                filter(self.ASSOC_par.mindetpop==self.mindetpop).\
                filter(self.ASSOC_par.minarraypop==self.minarraypop).\
                filter(self.ASSOC_par.duration==self.duration).\
                one()
            self.passocid=self.Passoc_Q.passocid
            #embed()
        else:
            print('process already in table: Assoc params table')
            self.passocid=self.Passoc_Q.passocid
            print(self.Passoc_Q)

        self.db_connected=True

        try:
            self.Affiliation_Q=self.session.query(self.Affiliation).filter(self.Affiliation.net==self.net).all()
        except Exception as ex1:
            print('Error with network retrieving', ex1)
            exit(0)
        refSTA=[]
        for aai in self.Affiliation_Q:

            try:
                STA_dataM=self.session.query(self.Site).filter(self.Site.sta==aai.sta).one()

            except Exception as ex1:
                #print
                print('there is more than just one station:', aai.sta,'  ',ex1)
                embed()
                exit()
            #embed()
            refSTA.append(STA_dataM.refsta)

        refstations_l=list(set(refSTA))
        refsta=[]
        #embed()
        for aai in refstations_l:
            STA_dataM=self.session.query(self.Site).filter(self.Site.refsta==str(aai)).all()
            array_lo=[]
            array_la=[]
            array_el=[]
            for sta_i in STA_dataM:
                array_la.append(sta_i.lat)
                array_lo.append(sta_i.lon)
                array_el.append(sta_i.elev)
            array_la=np.asarray(array_la)
            array_lo=np.asarray(array_lo)
            array_el=np.asarray(array_el)
            refsta.append({'lon':np.mean(array_lo),'lat':np.mean(array_la),'elev':np.mean(array_el),'name':aai,'numsta':len(array_la)})
        self.det_tot=[]
        self.fdtable_name=[]
        self.fktables_names=[]
        for aai in refsta:
            #print('getting data from:',aai['name'])
            # here it looks for any detections in all the tables
            for ti in range(self.num_tables):
                try:

                    fd_res=self.session.query(self.fdtables[ti]).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).all()
                    times_ini=self.session.query(self.fdtables[ti].timeini).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).all()
                    fk_table_names=self.session.query(self.fdtables[ti].fktablename).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).all()
                    #embed()
                    if len(fk_table_names)>0:
                        for tnfk in fk_table_names:
                            self.fktables_names.append(tnfk[0])

                            #embed()
                            '''
                            try:
                                if (tnfk[0] in self.dict_namefk)==False:
                                            self.fdtables=[]

                                    self.dict_namefk[tnfk[0]]=type(str(tnfk[0]) ,(schema.fk_results,),{'__tablename__':str(tnfk[0])})

                            except Exception as ex1:
                                print(ex1,'303')
                            '''

                except Exception as x1:
                    print('There is an error',x1)
                    embed()
                    exit()

        tt = np.unique(self.fktables_names)
        self.fktables=[]
        for  fki in tt:
            self.fktables.append(type(str(fki),(schema.fk_results,),{'__tablename__':str(fki)}))

        return self.db_connected

    def data_retrievingS(self,win_start, win_end):

        self.time_initial=win_start
        self.time_end=win_end

        id_res=self.session.query(self.ASSOC_results).filter(self.ASSOC_results.net==self.net).filter(self.ASSOC_results.passocid ==self.passocid).filter(self.ASSOC_results.timeini>self.time_initial).all()
        if len(id_res)>1:
            return 1
        else:
            print('start getting data for analysis')
        self.Detection=[]
        self.Detection_Q=[]
        #embed()
        refSTA=[]
        for aai in self.Affiliation_Q:

            try:
                STA_dataM=self.session.query(self.Site).filter(self.Site.sta==aai.sta).one()

            except Exception as ex1:
                #print
                print('there is more than just one station:', aai.sta,'  ',ex1)
                embed()
                exit()
            #embed()
            refSTA.append(STA_dataM.refsta)

        refstations_l=list(set(refSTA))
        refsta=[]
        #embed()
        for aai in refstations_l:
            STA_dataM=self.session.query(self.Site).filter(self.Site.refsta==str(aai)).all()
            array_lo=[]
            array_la=[]
            array_el=[]
            for sta_i in STA_dataM:
                array_la.append(sta_i.lat)
                array_lo.append(sta_i.lon)
                array_el.append(sta_i.elev)
            array_la=np.asarray(array_la)
            array_lo=np.asarray(array_lo)
            array_el=np.asarray(array_el)
            refsta.append({'lon':np.mean(array_lo),'lat':np.mean(array_la),'elev':np.mean(array_el),'name':aai,'numsta':len(array_la)})
        self.det_tot=[]
        self.fdtable_name=[]
        '''
        self.fktables_names=[]
        for aai in refsta:
            #print('getting data from:',aai['name'])
            # here it looks for any detections in all the tables
            for ti in range(self.num_tables):
                try:

                    fd_res=self.session.query(self.fdtables[ti]).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).filter(self.fdtables[ti].timeini>=self.time_initial).filter(self.fdtables[ti].timeini<=self.time_end).all()
                    times_ini=self.session.query(self.fdtables[ti].timeini).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).filter(self.fdtables[ti].timeini>=self.time_initial).all()
                    fk_table_names=self.session.query(self.fdtables[ti].fktablename).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).filter(self.fdtables[ti].timeini>=self.time_initial).filter(self.fdtables[ti].timeini<=self.time_end).distinct().all()
                    #embed()
                    if len(fk_table_names)>0:
                        for tnfk in fk_table_names:
                            self.fktables_names.append(tnfk[0])


                except Exception as x1:
                    print('There is an error',x1)
                    embed()
                    exit()

        self.fktables=[]
        for  fki in self.fktables_names:
            self.fktables.append(type(str(fki),(schema.fk_results,),{'__tablename__':str(fki)}))

        '''
        for aai in refsta:
            print('getting data from:',aai['name'])
            for ti in range(self.num_tables):
                try:

                    fd_res=self.session.query(self.fdtables[ti]).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).filter(self.fdtables[ti].timeini>=self.time_initial).filter(self.fdtables[ti].timeini<=self.time_end).all()
                    times_ini=self.session.query(self.fdtables[ti].timeini).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).filter(self.fdtables[ti].timeini>=self.time_initial).all()
                    if len(fd_res)>0:
                        print('length results:', len(fd_res))

                        times_ini=np.asarray(times_ini)
                        for tt in range(self.num_tables):
                        #fktable=self.dict_namefk[fk_table_names[0][ti]]
                            Query_fktempMA_all=self.session.query(self.fktables[tt]).filter(self.fktables[tt].pfkid==self.pfkid).filter(self.fktables[tt].timeini>=np.min(times_ini)).filter(self.fktables[tt].timeini<=np.max(times_ini)).filter(self.fktables[tt].sta==aai['name']).all()
                        #Query_fktempMA_all_timeini=self.session.query(fktable.timeini).filter(fktable.pfkid==self.pfkid).filter(fktable.timeini>=np.min(times_ini)).filter(fktable.timeini<=np.max(times_ini)).filter(fktable.sta==aai['name']).all()
                            Query_fktempMA_all_timeini=self.session.query(self.fktables[tt].timeini).filter(self.fktables[tt].pfkid==self.pfkid).filter(self.fktables[tt].timeini>=self.time_initial).filter(self.fktables[tt].timeini<=self.time_end).filter(self.fktables[tt].sta==aai['name']).all()
                            qt_all=[]
                            for qt in Query_fktempMA_all_timeini:
                                qt_all.append(qt[0])
                            qt_all=np.asarray(qt_all)
                            for dqi in range(len(fd_res)):
                            #for dqi in fd_res:
                                try:
                                    res_ind=np.where(qt_all==fd_res[dqi].timeini)
                                    Query_fktempMAX=Query_fktempMA_all[res_ind[0][0]]
                                    x=True
                                except Exception as x1:
                                    x=False
                                if x == True:
                            #det1 = inf_det_global(aai['lat'], aai['lon'], UTCDateTime(Query_fktempMAX.timeini).datetime, Query_fktempMAX.bz,  Query_fktempMAX.fval, aai['numsta'],DetID=dqi.fdid)
                                    det1 = (aai['lat'], aai['lon'], fd_res[dqi].timeini, Query_fktempMAX.bz,  Query_fktempMAX.fval, aai['numsta'], fd_res[dqi].fdid,aai['name'],self.fdtables_names[ti])
                                    self.det_tot.append(det1)

                                else:
                                    continue

                except Exception as x1:
                    print('There is an error',x1)
                    embed()
                    exit()
        print('all data retrieved')
        return 2


    def data_processingASSOC(self,t_start,t_end,src_win,max_prop_tm):
        '''

        '''
        print('data processing',short_time(UTCDateTime(self.time_initial)),short_time(UTCDateTime(self.time_end)))
        det_list = data_io.db2dets(self.det_tot)
        min_array_pop=self.minarraypop
        EVIDs=[]
        if len(det_list)>1:
            try:
                events = []
                event_qls = []
                window_start=[]
                window_end=[]
                duration_dd = int((t_end - t_start).astype('m8[s]').astype(float) / 60.0)
                #duration_dd = int((t_end - t_start) / 60.0)
                for dt in range(0, duration_dd, int(src_win)):
                    win_start = t_start + np.timedelta64(dt, 'm')
                    win_end = t_start + np.timedelta64(dt + int(src_win + max_prop_tm), 'm')
                    print('\n' + "Computing associations for:", win_start, " - ", win_end)
                    temp = [(n, det) for n, det in enumerate(det_list) if np.logical_and(win_start <= det.peakF_UTCtime, det.peakF_UTCtime <= win_end)]
                    key = [pair[0] for pair in temp]
                    new_list = [pair[1] for pair in temp]

                    # run analysis
                    if len(new_list)>1:
                        if self.trimthresh=='None':
                            self.trimthresh=None
                        labels, dists = hjl.run(new_list, self.clusterthresh, dist_max=self.distmax, bm_width=self.beamwidth, rng_max=self.rangemax,  pool=self.pl,trimming_thresh=self.trimthresh)
                        clusters, qualities = hjl.summarize_clusters(labels, dists,population_min=int(self.mindetpop))
                        for n in range(len(clusters)):
                            events += [[key[n] for n in clusters[n]]]
                            event_qls += [10.0**(-qualities[n])]
                            window_start.append(UTCDateTime(win_start.astype(datetime)).timestamp)
                            window_end.append(UTCDateTime(win_end.astype(datetime)).timestamp)
                event_cnt = len(events)
                for n1 in range(event_cnt):
                    for n2 in range(n1 + 1, event_cnt):
                        if len(events[n1]) > 0 and len(events[n2]) > 0:
                            set1, set2 = set(events[n1]), set(events[n2])
                            rel_overlap = len(set1.intersection(set2)) / min(len(set1), len(set2))

                            if rel_overlap > 0.5:
                                events[n1], events[n2] = list(set1.union(set2)), []
                                event_qls[n1], event_qls[n2] = max(event_qls[n1], event_qls[n2]), -1.0
                for n, ev_ids in enumerate(events):
                    if len(ev_ids) > 0:
                        locs = np.array([[det_list[j].latitude, det_list[j].longitude] for j in ev_ids])
                        #embed()
                        unique_cnt = max(len(np.unique(locs[:, 0])), len(np.unique(locs[:, 1])))
                        if unique_cnt < int(min_array_pop):
                            events[n] = []
                            event_qls[n] = -1.0
                events = [ei for ei in events if len(ei) > 0]
                event_qls = [eqi for eqi in event_qls if eqi > 0]
                print("Identified events and qualities:")
                for n in range(len(events)):
                    print('\t', events[n], '\t', event_qls[n])
                lastEVENTIDQ=self.session.query(func.max(self.ASSOC_results.eventid)).all()
                lastEVENTID=lastEVENTIDQ[0][0]
                if lastEVENTID is None:
                    lastEVENTID=int(0)

                for n in range(len(events)):
                    for nn in range(len(events[n])):
                        det_id = events[n][nn]
                        id_res=self.session.query(self.ASSOC_results).filter(self.ASSOC_results.net==self.net)\
                                                  .filter(self.ASSOC_results.fdid ==self.det_tot[det_id][6])\
                                                  .filter(self.ASSOC_results.passocid ==self.passocid)\
                                                  .filter(self.ASSOC_results.timeini==self.time_initial)\
                                                  .filter(self.ASSOC_results.timeend==self.time_end)\
                                                  .filter(self.ASSOC_results.qdetcluster==event_qls[n])\
                                                  .filter(self.ASSOC_results.fdtable==self.det_tot[det_id][8])\
                                                  .filter(self.ASSOC_results.sta==self.det_tot[det_id][7]).all()
                        id_resC=self.session.query(self.ASSOC_results).count()+1

                        if bool(id_res)==False:

                            res=self.ASSOC_results(associd=id_resC,\
                                    fdid=self.det_tot[det_id][6],\
                                    eventid=int(lastEVENTID+1+n),\
                                    passocid=self.passocid,\
                                    net=self.net,\
                                    timeini=window_start[n],\
                                    timeend=window_end[n],\
                                    qdetcluster=event_qls[n],\
                                    fdtable=self.det_tot[det_id][8],\
                                    sta=self.det_tot[det_id][7])
                            self.session.add(res)
                            self.session.commit()
                print('associations written', len(events))

            except Exception as ex1:
                print('error running assoc:',ex1)
                embed()
                exit()

            '''

            for ii in range(len(EVIDsS)):
                EVIDsN.append({'ID':self.det_tot[sorted_index[ii]].ID,'eventID':EVIDsS[ii],'qdetcluster':DAQS[ii],'qassoc':CAQS[ii],'fdtable':self.fdtable_nameS[ii]})
            lastEVENTIDQ=self.session.query(func.max(self.ASSOC_results.eventid)).all()
            lastEVENTID=lastEVENTIDQ[0][0]
            if lastEVENTID is None:
                lastEVENTID=int(0)

            for ev1 in EVIDsN:
                id_res=self.session.query(self.ASSOC_results).filter(self.ASSOC_results.net==self.net)\
                                          .filter(self.ASSOC_results.fdid ==ev1['ID'])\
                                          .filter(self.ASSOC_results.passocid ==self.passocid)\
                                          .filter(self.ASSOC_results.timeini==self.time_initial)\
                                          .filter(self.ASSOC_results.timeend==self.time_end)\
                                          .filter(self.ASSOC_results.qdetcluster==ev1['qdetcluster'])\
                                          .filter(self.ASSOC_results.fdtable==ev1['fdtable'])\
                                          .filter(self.ASSOC_results.qassoc==ev1['qassoc']).all()
                id_resC=self.session.query(self.ASSOC_results).count()+1

                if bool(id_res)==False:

                    res=self.ASSOC_results(associd=id_resC,\
                            fdid=ev1['ID'],\
                            eventid=int(lastEVENTID+ev1['eventID']),\
                            passocid=self.passocid,\
                            net=self.net,\
                            timeini=self.time_initial,\
                            timeend=self.time_end,\
                            qdetcluster=ev1['qdetcluster'],\
                            fdtable=ev1['fdtable'],\
                            qassoc=ev1['qassoc'])
                    self.session.add(res)
                    self.session.commit()
        print('associations written', len(EVIDs))


        id_resC=self.session.query(self.ASSOC_results).count()+1
        #embed()
        lastEVENTIDQ=self.session.query(func.max(self.ASSOC_results.eventid)).all()
        lastEVENTID=(lastEVENTIDQ[0][0])

        if lastEVENTID is None:
            lastEVENTID=int(0)

        res=self.ASSOC_results(associd=id_resC,\
                passocid=self.passocid,\
                net=self.net,\
                fdtable=self.net,\
                timeini=self.time_initial,\
                timeend=self.time_end,\
                #eventid=lastEVENTID+bytes(1))
                eventid=lastEVENTID)
        self.session.add(res)
        self.session.commit()

'''

    def data_processing(self):
        '''
        Constructor
        '''
        try:
            self.Affiliation_Q=self.session.query(self.Affiliation).filter(self.Affiliation.net==self.net).all()
        except Exception as ex1:
            print('Error with network retrieving', ex1)
            exit(0)

        for day in np.arange(self.dayofyearini,self.dayofyearend):

            print('year:',self.year,'day:', day)
            self.dayofyear=day
            t_start = UTCDateTime(year=self.year, julday = day, hour=0, minute=0)
            t_start = np.datetime64(t_start,'s')
            t_end = UTCDateTime(year=self.year, julday = day, hour=23, minute=59)
            t_end = np.datetime64(t_end,'s')
            days = self.jdayend-self.jdayini
            days = days * 24 * 60
            duration_dd = int((t_end - t_start).astype('m8[s]').astype(float) / 60.0)
            max_prop_tm = int((self.rangemax / 0.22) / self.duration)
            src_win = int(max_prop_tm * 0.5)
            win_start = UTCDateTime(t_start.astype(datetime)).timestamp
            win_end = UTCDateTime(t_end.astype(datetime)).timestamp
            retrieve_ret=self.data_retrievingS(win_start,win_end)
            for aa in range(len(self.det_tot)):
                print(UTCDateTime(self.det_tot[aa][2]))

            # retrieve time based on UTC date time string
            if retrieve_ret==2:
                self.data_processingASSOC(t_start,t_end,src_win,max_prop_tm)
            elif retrieve_ret==1 :
                print('time period already analyzed')
            else:
                print('there was an error or not data')


        #return self.data_processed


if __name__ == '__main__':
    pdetect=AssocInfraPy_LANL('../conf_files/InfraConfig_Assoc')
    pdetect.database_connecting()
    pdetect.data_processing()
