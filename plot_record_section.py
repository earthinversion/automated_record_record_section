import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os,yaml
from obspy import UTCDateTime, read
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
from obspy.core import Stream
from matplotlib.transforms import blended_transform_factory
import warnings
warnings.filterwarnings('ignore')

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
                # print()
                # print("event_description: ",dir(cat))
                # print("event_description: ",dir(cat.origins[0]))
                # print("event_description: ",cat.event_descriptions[0].text)
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
    # print(df.head())
    df.to_csv(catalogtxt, index=False)
    return df

## Reading yml file for input parameters
with open('input_file.yml') as f:
    inp = yaml.load(f, Loader=yaml.FullLoader)

minlat, maxlat = inp['minlat'], inp['maxlat']
minlon, maxlon = inp['minlon'], inp['maxlon']
starttime = inp['starttime']
endtime = inp["endtime"]
minmagnitude = inp["minmagnitude"]
maxmagnitude = inp["maxmagnitude"]
catalogxmlloc = inp["catalogxmlloc"]
client_name = inp["client_name"]
channel = inp["channel"]
location = inp["location"]
record_section_relative_beg_time = inp["record_section_relative_beg_time"]
record_section_relative_end_time = inp["record_section_relative_end_time"]

# defining output catalog file location and name
catalogxml = catalogxmlloc+'events.xml'
catalogtxt = catalogxmlloc+'events.txt'
stationxml = 'stations.xml'
stationtxt = 'stations.txt'




# retrieve information for all available events during the given time range
kwargs = {'starttime': starttime, 'endtime': endtime, 
    'minlatitude': minlat, 'maxlatitude': maxlat, 
    'minlongitude': minlon, 'maxlongitude': maxlon,
    'minmagnitude': minmagnitude, 'maxmagnitude': maxmagnitude}

client = Client(client_name)
print("--> Obtaining catalog")
try:
    catalog = client.get_events(**kwargs)
except:
    print(f"ConnectionResetError while obtaining the events from the client - {client_name}")
catalog.write(catalogxml, 'QUAKEML') #writing xml catalog

df = extract_catalog_info_to_txt(catalog, catalogtxt) ## Writing catalog in txt file and returning pandas dataframe



### Download waveforms and add to the stream
eq_id = 0
eq_origin_time = UTCDateTime(df['eventTime'].values[eq_id])
ev_time_str = df['eventTime'].values[eq_id].split("T")[0]
wf_starttime = eq_origin_time - record_section_relative_beg_time
wf_endtime = eq_origin_time + record_section_relative_end_time
evLat = df['eventLat'].values[eq_id]
evLon = df['eventLon'].values[eq_id]

print("--> Obtaining stations information")
inventory = client.get_stations(startbefore=wf_starttime, endafter=wf_endtime,network="*", station="*", channel=channel, level='response', minlongitude=minlon, maxlongitude=maxlon, minlatitude=minlat, maxlatitude=maxlat, location=location)
inventory.write(stationtxt, 'STATIONTXT',level='station')

df_stationlist = pd.read_csv(stationtxt,sep="|")
# df_stationlist = df_stationlist.iloc[16:100,:]
dists = []
print("--> Downloading waveforms for the event: {} @{}".format(df['eventDescrp'].values[eq_id],eq_origin_time))
st = Stream()
count = 0
for ii, net, sta, lat, lon in zip(np.arange(df_stationlist.shape[0]), df_stationlist['#Network'],df_stationlist['Station'],df_stationlist['Latitude'],df_stationlist['Longitude']):
    try:
        tr = client.get_waveforms(net, sta, location, channel, wf_starttime, wf_endtime, attach_response = True)
        dist, baz, _ = gps2dist_azimuth(lat,lon,evLat,evLon)
        

        st += tr
        dists.append(dist)
        st[count].stats['distance'] = dist

        print(f"{ii+1}/{df_stationlist.shape[0]} Successfully downloaded trace for {net}-{sta}-{channel}, dist: {dist/1000:.2f} km")
        count += 1
    except:
        print(f"{ii}/{df_stationlist.shape[0]} Failed {net}-{sta}-{channel}: ({lat:.2f}, {lon:.2f})")

    #plot N number of traces only to increase visibility
    if len(st)>=inp['max_num_trace_to_plot']:
        break

# bandpass filter
if inp['filter_trace']:
    print("--> Bandpass filtering")
    st.filter('bandpass', freqmin=inp['bandpass_filter']['freqmin'], freqmax=inp['bandpass_filter']['freqmax'])

# remove instrument response
st_rem = st.copy()
if inp['rem_response']:
    print("--> Removing response")
    st_rem.remove_response(output='DISP')

st_rem.write(f"stream_rs_{ev_time_str}.mseed", format="MSEED")

print("--> Plotting record section")

fig = plt.figure()
st_rem.plot(type='section',grid_linewidth=0.1,dist_degree=False, fig=fig, time_down=True)

ax = fig.axes[0]
transform = blended_transform_factory(ax.transData, ax.transAxes)
for tr in st_rem:
    ax.text(tr.stats.distance / 1e3, 1.0, tr.stats.station, rotation=270,
            va="bottom", ha="center", transform=transform, zorder=10)
plt.title(f"Record section for {df['eventMag'].values[eq_id]} {df['magUnit'].values[eq_id]} EQ @ {df['eventDescrp'].values[eq_id]}\n\n", fontsize=14)
plt.savefig(f"record_section_{ev_time_str}.png",dpi=300,bbox_inches='tight')