from astropy.table import Table
import numpy as np
from kadi import events
from Ska.engarchive import fetch
from mica.archive import aca_l0
from Chandra.Time import DateTime

# only look for fids in slots 0, 1, 2
FID_SLOTS = [0, 1, 2]

COLS = []
SLOT_COLS = ['AOIMAGE', 'AOACFCT']
ALL_COLS = COLS + [field + '%s' % slot
                   for field in SLOT_COLS for slot in FID_SLOTS]

DROP_PAD = 450
start = '1999:001'
dwells = events.dwells.filter(start=start)

print """Searching dwells from {} for fid drops which occur more
than {} secs before the end of the dwell.  Uses the first value
of AOIMAGE for the slot to determine if it should be a FID""".format(
    start, DROP_PAD)

# Make a dictionary to store fid drops; use dwell start time as key
all_drop = {}
for dwell in dwells:

    if (dwell.start == '2015:076:03:53:33.193'):
        print """
SCS107 for 16344, which isn't caught by this script
because the fids are all lost within the last {} secs
of the observation""".format(DROP_PAD)
        continue

    actual_fid_slots = []
    pcad_data = fetch.MSIDset(ALL_COLS, dwell.tstart, dwell.stop, filter_bad=True)
    # Use the first value of AOIMAGE to determine if the slot was ever a FID
    for fid_slot in FID_SLOTS:
        if pcad_data['AOIMAGE{}'.format(fid_slot)].vals[0] == 'FID ':
            actual_fid_slots.append(fid_slot)

    # For each fid slot, check to see if we're still tracking the fid up to
    # the end of the dwell.  Record loss times if fid lost more than DROP_PAD
    # before the end of the dwell.
    drop_slot_time = {}
    drop_slot_dtime = {}
    for slot in actual_fid_slots:
        if pcad_data['AOIMAGE{}'.format(slot)].vals[-1] == 'FID ':
            continue
        if pcad_data['AOACFCT{}'.format(slot)].vals[-1] == 'TRAK':
            continue
        if np.any(pcad_data['AOACFCT{}'.format(slot)].vals != 'TRAK'):
            time = pcad_data['AOACFCT{}'.format(slot)].times[
                np.flatnonzero(
                    pcad_data['AOACFCT{}'.format(slot)].vals != 'TRAK')[0]]
            # Some HRC observations lose all fids ~400 seconds before the end of the nominal dwell
            # We aren't explicitly interested in these cases for this review.
            if time < (DateTime(dwell.stop).secs - DROP_PAD):
                drop_slot_dtime[slot] = "{:2f}".format(DateTime(dwell.stop).secs - time)
                drop_slot_time[slot] = "%.2f" % time
        else:
            raise ValueError

    # If we lost a fid, print out some details
    if len(drop_slot_time):
        obsid = None
        try:
            obsid = dwell.get_obsid()
            print "\nFids lost for dwell at {} Obsid {}".format(dwell.start, obsid)
        except:
            print "\nFids lost for dwell at {}".format(dwell.start)
            pass
        for slot in drop_slot_time:
            print "Fid {} lost at {} (dwell stop - {} secs)".format(slot,
                                                                    drop_slot_time[slot],
                                                                    drop_slot_dtime[slot])
        # if more than one fid lost, this was "interesting" and we should know why
        if len(drop_slot_time) > 1:
            # HRC shutter test 1201
            if (dwell.start == '1999:242:23:31:59.988'):
                print "HRC shutter test 1201"
                continue
            # HRC shutter test 1211 (processed with 1 fid)
            if (dwell.start == '1999:245:05:13:05.795'):
                print "HRC shutter test 1211 processed 1 fid"
                continue
            # HRC shutter test 1212 (processed with no fids)
            if (dwell.start == '1999:245:05:58:40.495'):
                print "HRC shutter test 1212 processed no fids"
                continue
            # HRC shutter test 1214 (processed with 1 fid)
            if (dwell.start == '1999:245:07:32:01.095'):
                print "HRC shutter test 1214 processed 1 fid"
                continue
            # HRC shutter test 1216
            if (dwell.start == '1999:245:09:05:25.795'):
                print "HRC shutter test 1216"
                continue
            # Obsid 1217 Two fids to start, one not used
            if (dwell.start == '1999:245:09:51:57.895'):
                print "Obsid 1217 only two fids to start, one not used"
                continue
            if (dwell.start == '1999:263:11:34:23.551'):
                print "Obsid 361 fid drop slot 2"
                continue
            if (dwell.start == '1999:347:21:44:43.620'):
                print "Obsid 1507 SI Safing test"
                continue
            if dwell.start == '2013:076:07:42:19.121':
                print "{} all fids lost SCS107 (missing in kadi.events due to known bug)"
                continue
            if events.scs107s.filter(dwell.start, dwell.stop).count():
                print "{} all fids lost SCS107".format(obsid)
                continue
            if events.normal_suns.filter(dwell.start, dwell.stop).count():
                print "{} all fids lost NSM".format(obsid)
                continue
            if dwell.start == '2015:048:04:47:51.658':
                print "obsid 16091 two fid drop"
                continue
            # Toss an error if we have an unexplained multi-fid-drop
            raise ValueError

        all_drop[dwell.tstart] = drop_slot_time

