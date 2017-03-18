import ftplib
import tempfile
import tarfile
import os
import sys

def ymd_nexrad(filename):
    year = filename[4:8]
    month = filename[8:10]
    day = filename[10:12]
    hour = filename[13:15]
    minute = filename[15:17]
    second = filename[17:19]
    return year, month, day, hour, minute, second

def fetch_has_tar(ftp, filename, odir):
    fh = tempfile.NamedTemporaryFile()
    ftp.retrbinary('RETR ' + filename, fh.write)
    my_tar = tarfile.open(fh.name)
    try:
        fnames = my_tar.getnames()
        failed_state = False
    except tarfile.ReadError:
        failed_state = True
        inventory = [filename + ' FAILED']
    if not failed_state:
        inventory = []
        for fname in fnames:
            fhout = my_tar.extractfile(fname)
            year, month, day, hour, minute, second = ymd_nexrad(fname)
            pth = os.path.join(odir, year, month, day)

            try:
                os.makedirs(pth)
            except FileExistsError:
                pass #directory exists

            fhwrite = open(os.path.join(pth, fname), 'wb')
            fhwrite.writelines(fhout.readlines())
            fhwrite.close()
            fhout.close()
            inventory.append(pth + ',' + fname)
    return inventory

def fetch_has(has, odir, inv_dir=None, loud = False):
    if inv_dir is None:
        inv_dir = odir
    ftp = ftplib.FTP('ftp.ncdc.noaa.gov')
    ftp.login()
    invent = os.path.join(inv_dir, has+'_inventory.txt')
    invfh = open(invent, 'w')
    ftp.cwd('pub/has/'+has+'/')
    lst = ftp.nlst()
    for tarfile in lst:
        if loud:
            print('Doing', tarfile)
        this_inventory = fetch_has_tar(ftp, tarfile, odir)
        this_inventory_cr = [item + '\n' for item in this_inventory]
        invfh.writelines(this_inventory_cr)

    invfh.close()
    ftp.close()

if __name__ == '__main__':
    has = sys.argv[1]
    odir = sys.argv[2]
    invdir = None
    if len(sys.argv) == 4:
        invdir = sys.argv[3]

    fetch_has(has, odir, inv_dir=None, loud = False):



