import numpy as np


def get_mean_locations(session, net, Affiliation, Site, t0_E_jday, te_E_jday):
    """
    Get mean location information for a given network.

    Parameters
    ----------
    session : SQLAlchemy Session instance
    net : str
        Network code.
    Affiliation, Site : SQLAlchemy ORM mapped table classes
        Affiliation, Site classes (not instances).
    t0_E_jday, te_E_jday : int
        Julian date (YYYJJJ) Site.ondate, offdate.

    Returns
    -------
    refsta : list of dict
        One dict for each unique refsta, each dict contains mean lat, lon,
        elev, values in 'lat', 'lon', 'elev' keys, name of refsta in 'name'
        key, and number of stations in the means in 'numsta'.

    """
    refSTA = []

    try:
        Affiliation_Q = session.query(Affiliation).filter(Affiliation.net==net).all()
    except Exception as ex1:
        print('Error with network retrieving', ex1)
        exit(0)

    for aai in Affiliation_Q:
        try:
            q = session.query(Site).filter(Site.sta==aai.sta)
            q = q.filter(((Site.offdate==-1) | (Site.offdate>te_E_jday)) & (Site.ondate<t0_E_jday))
            STA_dataM = q.one()
        except Exception as ex1:
            print('there is more than just one station:', aai.sta,'  ',ex1)
            exit()
        refSTA.append(STA_dataM.refsta)

    refstations_l = list(set(refSTA))
    refsta = []
    for aai in refstations_l:
        q = session.query(Site).filter(Site.refsta==str(aai))
        q = q.filter(((Site.offdate>te_E_jday) | (Site.offdate==-1)) & (Site.ondate<t0_E_jday))
        STA_dataM = q.all()

        array_lo = []
        array_la = []
        array_el = []
        for sta_i in STA_dataM:
            array_la.append(sta_i.lat)
            array_lo.append(sta_i.lon)
            array_el.append(sta_i.elev)

        array_la = np.asarray(array_la)
        array_lo = np.asarray(array_lo)
        array_el = np.asarray(array_el)

        d = {'lon': np.mean(array_lo),
             'lat': np.mean(array_la),
             'elev': np.mean(array_el),
             'name': aai,
             'numsta': len(array_la)}
        refsta.append(d)

        return refsta
