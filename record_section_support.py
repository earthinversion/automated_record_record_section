import pandas as pd


def extract_catalog_info_to_txt(catalog, catalogtxt):
    '''
    Writing catalog in txt file and returning pandas dataframe
    :type catalog: :class: 'obspy.core.event.catalog.Catalog'
    :param catalog: earthquake catalog
    :type catalogtxt: string
    :param catalogtxt: output catalog text file name
    '''
    evtimes, evlats, evlons, evdps, evmgs, evmgtps, evdess = [], [], [], [], [], [], []
    for cat in catalog:
        try:
            try:
                evtime,evlat,evlon,evdp,evmg,evmgtp, evdes =cat.origins[0].time,cat.origins[0].latitude,cat.origins[0].longitude,cat.origins[0].depth/1000,cat.magnitudes[0].mag,cat.magnitudes[0].magnitude_type, cat.event_descriptions[0].text

            except:
                evtime,evlat,evlon,evdp,evmg,evmgtp, evdes=cat.origins[0].time,cat.origins[0].latitude,cat.origins[0].longitude,cat.origins[0].depth/1000,cat.magnitudes[0].mag,"-9999", cat.event_descriptions[0].text
            evtimes.append(str(evtime))
            evlats.append(float(evlat))
            evlons.append(float(evlon))
            evdps.append(float(evdp))
            evmgs.append(float(evmg))
            evmgtps.append(str(evmgtp))  
            evdess.append(evdes)      
        except:
            print(f"Unable to write for {evtime}")
    df = pd.DataFrame({'eventTime':evtimes, 'eventLat': evlats, 'eventLon': evlons, 'eventDepths': evdps, "eventMag" : evmgs, "magUnit": evmgtps, "eventDescrp": evdess})
    df = df.sort_values(by=['eventMag'], ascending=False)
    print(df.head())
    df.to_csv(catalogtxt, index=False)
    return df