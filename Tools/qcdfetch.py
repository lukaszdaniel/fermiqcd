#!/usr/bin/env python3
# -*- coding: utf-8 -*-

LOGO ="""
 #######   ######  ########     ######## ######## ########  ######  ##     ## 
##     ## ##    ## ##     ##    ##       ##          ##    ##    ## ##     ## 
##     ## ##       ##     ##    ##       ##          ##    ##       ##     ## 
##     ## ##       ##     ##    ######   ######      ##    ##       ######### 
##  ## ## ##       ##     ##    ##       ##          ##    ##       ##     ## 
##    ##  ##    ## ##     ##    ##       ##          ##    ##    ## ##     ## 
 ##### ##  ######  ########     ##       ########    ##     ######  ##     ## 
Created by Massimo Di Pierro - License GPL2 - all-to-all convertion utility
"""

USAGE ="""Usage:

  $ qcdfetch [options] source

Source can be a link to a NERSC ensable, a local file or a glob pattern 
(for example *.milc, folder/*.ildg)

Example: print help:
  $ qcdfetch -h

Example: download an ensable:
  $ qcdfetch source

Example: fetch an ensable in fermiqcd (.mdp) format.
  $ qcdfetch -c mdp source

Example: fetch an ensable in fermiqcd (.mdp) format an convert to float
  $ qcdfetch -c mdp -4 source

Example: fetch an ensable in fermiqcd (.mdp) format an convert to double
  $ qcdfetch -c mdp -8 source

Example: fetch an ensable in ildg (.ildg)
  $ qcdfetch -c ildg source

Example: fetch an ensable in fermiqcd format and break it into timeslices
  $ qcdfetch -c split.mdp source
  (will make one file per timeslice)

Example: fetch a propagator and convert in fermiqcd(.prop.mdp) format
  $ qcdfetch -c prop.mdp source

Example: fetch a propagator and convert in fermiqcd format and break into timeslices
  $ qcdfetch -c split.prop.mdp source
  (will make one file per timeslice)

Notice: qcdfetch runs in the folder where files are located unless, in case of 
downloads, a --destination folder is specified.

"""

##### imports #############################################################

from abc import ABC, abstractmethod
from ftplib import FTP
import urllib
import hashlib
import pickle
import os
import re
import sys
import time
from datetime import datetime
import optparse
import struct
import mmap
import glob
import io
import termios
import signal
import array
import fcntl
import logging
import traceback
import xml.dom.minidom as dom
import xml.parsers.expat as expat
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed

##### global variables #############################################################

NOW = datetime.now()
MAXBYTES = 1000  # Max number of bytes for buffered reading
PRECISION = {'f': 32, 'd': 64}
(X, Y, Z, T) = (1, 2, 3, 0)  # MDP index convention, used internally
DEFAULT_WIDGETS = [Percentage(), ' ', Bar()]

def notify(*args):
    """Wrapper for print function."""
    print(' '.join(map(str, args)))

# Logger setup for better message handling
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


##### class Lime #############################################################

class Lime:
    """
    LIME format reader and writer.
    Based on the LIME format specification: https://usqcd-software.github.io/c-lime/lime_1p2.pdf
    """

    @staticmethod
    def xml_parser(data):
        """
        Helper function to parse XML data.
        Args:
            data (str): XML content to parse.
        Returns:
            function: A function to retrieve XML tag content by name.
        """
        def f(name, dxml=dom.parseString(data)):
            return dxml.getElementsByTagName(name)[0].childNodes[0].nodeValue
        return f

    def __init__(self, filename, mode, version=1):
        """
        Initialize the Lime file reader/writer.

        Args:
            filename (str): Name of the LIME file.
            mode (str): Mode of file operation ('r' for read, 'w' for write).
            version (int, optional): Version of the LIME format. Defaults to 1.
        """
        self.magic = 1164413355
        self.version = version
        self.filename = filename
        self.mode = mode
        self.file = open(filename, mode + 'b')
        self.records = []

        if mode in ('r', 'rb'):
            self._read_records()

        self.dump_info()

    def _read_records(self):
        """
        Read LIME records from the file.
        """
        while True:
            header = self.file.read(144)
            if not header:
                break
            try:
                magic, _, size, name = struct.unpack('!iiq128s', header)
            except struct.error:
                raise IOError("Not in LIME format")

            if magic != self.magic:
                raise IOError("Not in LIME format")

            name = name.split(b'\0', 1)[0].decode()  # Clean up junk from file
            position = self.file.tell()
            self.records.append((name, position, size))  # In bytes

            padding = (8 - (size % 8)) % 8
            self.file.seek(size + padding, 1)  # Skip to next record

    def dump_info(self, filename=None):
        """
        Dump information about the records in the LIME file.

        Args:
            filename (str, optional): Name of the file to dump info to. Defaults to None.
        """
        with open(filename or (self.filename + '.fromlime.info'), 'w') as f:
            f.write("LIME records:\n")
            for name, pos, size in self.records:
                f.write(f'- {name} [{size} bytes]\n')
                if size < 1000:
                    self.file.seek(pos)
                    f.write('\n' + self.file.read(size).decode(errors='ignore') + '\n\n')

    def read(self, record):
        """
        Reads a LIME record by index.

        Args:
            record (int): Index of the record.

        Returns:
            tuple: (name, reader, size)
        """
        if self.mode not in ('r', 'rb'):
            raise RuntimeError("Reading is not supported in the current mode")

        name, position, size = self.records[record]
        self.file.seek(position)
        return name, self.file, size

    def __iter__(self):
        """
        Iterate through all LIME records.

        Yields:
            tuple: (name, reader, size)
        """
        for record in range(len(self.records)):
            yield self.read(record)

    def write(self, name, reader, size=None, chunk=MAXBYTES):
        """
        Write a LIME record.

        Args:
            name (str): Name of the record.
            reader (str or file-like object): Data to be written.
            size (int, optional): Size of the record. Defaults to None.
            chunk (int, optional): Maximum size of chunks to write at a time. Defaults to MAXBYTES.
        """
        if self.mode not in ('w', 'wb'):
            raise RuntimeError("Writing is not supported in the current mode")

        if isinstance(reader, str):
            if size is None:
                size = len(reader)
            reader = io.BytesIO(reader.encode())

        # Write record header
        position = self.file.tell()
        header = struct.pack('!iiq128s', self.magic, self.version, size, name.encode())
        self.file.write(header)

        # Write data in chunks
        self._write_data(reader, size, chunk)

        # Add padding bytes
        padding = (8 - (size % 8)) % 8
        self.file.write(b'\0' * padding)

        self.records.append((name, size, position))

    def _write_data(self, reader, size, chunk):
        """
        Helper function to write data from the reader to the file in chunks.

        Args:
            reader (file-like object): Data source.
            size (int): Size of the data.
            chunk (int): Chunk size for writing.
        """
        if hasattr(reader, 'read'):
            for _ in range(size // chunk):
                data = reader.read(chunk)
                if len(data) != chunk:
                    raise IOError("Unexpected EOF")
                self.file.write(data)

            # Write remaining data
            remaining = size % chunk
            data = reader.read(remaining)
            if len(data) != remaining:
                raise IOError("Unexpected EOF")
            self.file.write(data)
        else:
            for data in reader:
                self.file.write(data.encode() if isinstance(data, str) else data)

    def close(self):
        """Close the file."""
        self.file.close()

    def __len__(self):
        """
        Return the number of LIME records.

        Returns:
            int: Number of records.
        """
        return len(self.records)

    def keys(self):
        """
        Return the names of all LIME records.

        Returns:
            list: List of record names.
        """
        return [name for name, position, size in self.records]

def test_lime():
    """Test creating and reading a LIME file."""
    notify("Making a dummy LIME file and writing junk in it...")
    lime = Lime('test_lime_0.lime', 'w')
    lime.write('record1', '01234567')
    lime.write('record2', '012345678')
    file = io.BytesIO(b'0123456789')  # In-memory file
    lime.write('record3', file, 10)  # Write file content as record
    lime.close()

    notify("Reading the file back...")
    lime = Lime('test_lime_0.lime', 'r')
    notify(f"File contains {len(lime)} records")
    notify(f"They have names: {lime.keys()}")

    for name, reader, size in lime:
        notify(f"Record name: {name}\nRecord size: {size}\nRecord data: {reader.read(size).decode(errors='ignore')}")
    notify("LIME testing finished")
    lime.close()
    os.system('rm test_lime_0.*')


##### reunitarize #############################################################

def reunitarize(s):
    """
    Reunitarizes a list of complex numbers.

    Args:
        s: List of 12 complex numbers (a1re, a1im, ..., b3im).

    Returns:
        A new list with 18 reunitarized complex numbers.
    """
    (a1re, a1im, a2re, a2im, a3re, a3im, b1re, b1im, b2re, b2im, b3re, b3im) = s

    c1re = a2re * b3re - a2im * b3im - a3re * b2re + a3im * b2im
    c1im = -(a2re * b3im + a2im * b3re - a3re * b2im - a3im * b2re)
    c2re = a3re * b1re - a3im * b1im - a1re * b3re + a1im * b3im
    c2im = -(a3re * b1im + a3im * b1re - a1re * b3im - a1im * b3re)
    c3re = a1re * b2re - a1im * b2im - a2re * b1re + a2im * b1im
    c3im = -(a1re * b2im + a1im * b2re - a2re * b1im - a2im * b1re)

    return (a1re, a1im, a2re, a2im, a3re, a3im,
            b1re, b1im, b2re, b2im, b3re, b3im,
            c1re, c1im, c2re, c2im, c3re, c3im)

def check_unitarity(items, tolerance=1.0):
    """
    Checks if all numbers in the list are within the range [-tolerance, +tolerance].

    Args:
        items: List of values to check.
        tolerance: The acceptable range for values.

    Raises:
        RuntimeError: If any value is outside the acceptable range.
    """
    if any(x < -tolerance or x > tolerance for x in items):
        raise RuntimeError("Matrix is not unitary")

def reorder(data, order1, order2):
    """
    Reorders a list of complex numbers based on given orderings.

    Args:
        data: List of complex numbers.
        order1: The first ordering of indices.
        order2: The second ordering of indices.

    Returns:
        A reordered list of data.
    """
    k = len(data)  # 4*9*2
    m = len(order1)  # 4
    n = k // m  # 9*2

    items = [None] * k
    for i in range(k):
        items[n * order1[i // n] + i % n] = data[i]

    items = [items[n * order2[i // n] + i % n] for i in range(k)]
    return items

assert ''.join(reorder('AABBCCDD', [X, Y, Z, T], [Z, Y, X, T] )) == 'CCBBAADD'

##### Field readers #############################################################

class QCDFormat:
    site_order = [T, Z, Y, X]  # Always unused but for reference
    link_order = [X, Y, Z, T]  # Order of links at the site level

    def __init__(self, precision='f', endianess='>', is_gauge=True):
        """
        Constructor to initialize defaults.

        Args:
            precision: Precision ('f' or 'd').
            endianess: Endianess of the data.
            is_gauge: Whether the data is gauge data.
        """
        self.precision = precision
        self.endianess = endianess
        self.is_gauge = is_gauge

    def unpack(self, data):
        """
        Unpacks a string of bytes into a list of float/double numbers.

        Args:
            data: The data to unpack.

        Returns:
            A list of unpacked numbers.
        """
        if self.precision.lower() == 'f':
            n = len(data) // 4
        elif self.precision.lower() == 'd':
            n = len(data) // 8
        else:
            raise IOError("Incorrect input precision")

        items = struct.unpack(self.endianess + str(int(n)) + self.precision, data)

        if self.is_gauge:
            items = reorder(items, self.link_order, [T, X, Y, Z])
            check_unitarity(items)

        return items

    def pack(self, items):
        """
        Packs a list of numbers into a string of bytes.

        Args:
            items: The list of numbers to pack.

        Returns:
            A packed byte string.
        """
        if self.is_gauge:
            items = reorder(items, [T, X, Y, Z], self.link_order)

        n = len(items)
        return struct.pack(self.endianess + str(n) + self.precision, *items)

    @abstractmethod
    def read_header(self):
        """Reads file header or fails."""
        pass

    @abstractmethod
    def read_data(self, t, x, y, z):
        """Random access read."""
        pass

    @abstractmethod
    def write_header(self, precision, nt, nx, ny, nz):
        """Write header for new file."""
        pass

    @abstractmethod
    def write_data(self, data):
        """Write next site variables, in order."""
        pass

    def close(self):
        """Closes the file."""
        self.file.close()

class GaugeCold(QCDFormat):
    def __init__(self, nt=8, nx=4, ny=4, nz=4, precision='f', endianess='>', is_gauge=True):
        super().__init__(precision, endianess, is_gauge)
        self.lattice_size = (nt, nx, ny, nz)

    def read_header(self):
        nt, nx, ny, nz = self.lattice_size
        return (self.precision, nt, nx, ny, nz)

    def read_data(self, t, x, y, z):
        return [1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

class GaugeMDP(QCDFormat):
    site_order = [T, X, Y, Z]
    link_order = [T, X, Y, Z]

    def __init__(self, filename):
        super().__init__(None, '<', True)
        self.filename = filename
        self.header_format = '<60s60s60sLi10iii'
        self.header_size = 60 + 60 + 60 + 14 * 4  # Sum of sizes in header_format
        self.offset = None
        self.site_size = None
        self.base_size = 4 * 9 * 2  # Base size for site (float precision)
        self.lattice_size = None

    def read_header(self):
        """
        Reads the header of the data file, determining precision and lattice size.

        Returns:
            Tuple: The precision, nt, nx, ny, nz of the lattice.
        """
        with open(self.filename, 'rb') as f:
            header = f.read(self.header_size)
            items = struct.unpack(self.header_format, header)

            if items[3] != 1325884739:
                notify("Warning: This does not appear to be an MDP file, but it could be wrong.")

            nt, nx, ny, nz = items[5:9]
            self.site_size = items[15]

            if self.site_size == self.base_size * 4:
                self.precision = 'f'
            elif self.site_size == self.base_size * 8:
                self.precision = 'd'
            else:
                raise IOError("Unable to determine input precision")

            self.offset = f.tell()
            self.lattice_size = (nt, nx, ny, nz)

        return (self.precision, nt, nx, ny, nz)

    def write_header(self, precision, nt, nx, ny, nz):
        """
        Writes the header to the file in the LIME format.

        Args:
            precision: Precision ('f' or 'd').
            nt, nx, ny, nz: Lattice dimensions.
        """
        with open(self.filename, 'wb') as f:
            self.site_size = self.base_size * (4 if precision == 'f' else 8)
            now = NOW.isoformat()
            data = struct.pack(self.header_format.encode('utf-8'),
                               'File Type: MDP FIELD'.encode('utf-8'),
                               'none'.encode('utf-8'), #self.filename.encode('utf-8'),
                               now.encode('utf-8'),
                               1325884739, 4, nt, nx, ny, nz, 0, 0, 0, 0, 0, 0,
                               self.site_size, nt * nx * ny * nz)
            f.write(data)
            self.lattice_size = (nt, nx, ny, nz)
            self.precision = precision
            self.offset = f.tell()

    def read_data(self, t, x, y, z):
        """
        Reads the gauge data for a specific lattice point.

        Args:
            t, x, y, z: Coordinates of the lattice site.

        Returns:
            The unpacked data for the given coordinates.
        """
        nt, nx, ny, nz = self.lattice_size
        index = self.offset + (z + nz * (y + ny * (x + nx * t))) * self.site_size

        with open(self.filename, 'rb') as f:
            f.seek(index)
            data = f.read(self.site_size)

        return self.unpack(data)

    def write_data(self, data, target_precision=None):
        """
        Writes the given data to the file.

        Args:
            data: The data to be written.
            target_precision: Precision of the data (optional).
        """
        if len(data) != self.base_size:
            raise RuntimeError("Invalid data size")

        with open(self.filename, 'ab') as f:  # Open file in append mode
            f.write(self.pack(data))

    def convert_from(self, other, target_precision=None):
        """
        Converts data from another object.

        Args:
            other: Another object instance.
            target_precision: Desired precision for the conversion (optional).
        """
        precision, nt, nx, ny, nz = other.read_header()
        print(f"  (precision: {precision}, size: {nt}x{nx}x{ny}x{nz})")

        self.write_header(target_precision or precision, nt, nx, ny, nz)

        pbar = ProgressBar(widgets=DEFAULT_WIDGETS, maxval=nt).start()

        for t in range(nt):
            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        data = other.read_data(t, x, y, z)
                        self.write_data(data)

            pbar.update(t)

        pbar.finish()

class GaugeMDPSplit(GaugeMDP):
    def convert_from(self, other, target_precision=None):
        precision, nt, nx, ny, nz = other.read_header()
        print(f"  (precision: {precision}, size: {nt}x{nx}x{ny}x{nz})")

        pbar = ProgressBar(widgets=DEFAULT_WIDGETS, maxval=nt).start()

        for t in range(nt):
            slice_filename = self.filename.replace('split.prop.mdp', f't{t:04d}.prop.mdp')
            slice = GaugeMDP(slice_filename)
            slice.write_header(target_precision or precision, 1, nx, ny, nz)

            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        data = other.read_data(t, x, y, z)
                        slice.write_data(data)

            slice.close()
            pbar.update(t)

        pbar.finish()

class PropagatorMDP(QCDFormat):
    site_order = [T, X, Y, Z]
    link_order = [T, X, Y, Z]

    def __init__(self, filename):
        super().__init__(None, '<', False)
        self.filename = filename
        self.header_format = '<60s60s60sLi10iii'
        self.header_size = 60 + 60 + 60 + 14 * 4  # Sum of sizes in header_format
        self.offset = None
        self.site_size = None
        self.base_size = 16 * 9 * 2
        self.lattice_size = None

    def read_header(self):
        with open(self.filename, 'rb') as f:
            header = f.read(self.header_size)
            items = struct.unpack(self.header_format, header)

            if items[3] != 1325884739:
                notify("Warning: This does not appear to be an MDP file, but it could be wrong.")

            nt, nx, ny, nz = items[5:9]
            self.site_size = items[15]

            if self.site_size == self.base_size * 4:
                self.precision = 'f'
            elif self.site_size == self.base_size * 8:
                self.precision = 'd'
            else:
                raise IOError("Unable to determine input precision")

            self.offset = f.tell()
            self.lattice_size = (nt, nx, ny, nz)

        return (self.precision, nt, nx, ny, nz)

    def write_header(self, precision, nt, nx, ny, nz):
        with open(self.filename, 'wb') as f:
            self.site_size = self.base_size * (4 if precision == 'f' else 8)
            now = NOW.isoformat()
            data = struct.pack(self.header_format.encode('utf-8'),
                               'File Type: MDP FIELD'.encode('utf-8'),
                               'none'.encode('utf-8'), #self.filename.encode('utf-8'),
                               now.encode('utf-8'),
                               1325884739, 4, nt, nx, ny, nz, 0, 0, 0, 0, 0, 0,
                               self.site_size, nt * nx * ny * nz)
            f.write(data)
            self.lattice_size = (nt, nx, ny, nz)
            self.precision = precision
            self.offset = f.tell()

    def read_data(self, t, x, y, z):
        nt, nx, ny, nz = self.lattice_size
        index = self.offset + (z + nz * (y + ny * (x + nx * t))) * self.site_size

        with open(self.filename, 'rb') as f:
            f.seek(index)
            data = f.read(self.site_size)

        return self.unpack(data)

    def write_data(self, data, target_precision=None):
        if len(data) != self.base_size:
            raise RuntimeError("Invalid data size")

        with open(self.filename, 'ab') as f:  # Open file in append mode
            f.write(self.pack(data))

    def convert_from(self, other, target_precision=None):
        precision, nt, nx, ny, nz = other.read_header()
        print(f"  (precision: {precision}, size: {nt}x{nx}x{ny}x{nz})")
        self.write_header(target_precision or precision, nt, nx, ny, nz)

        pbar = ProgressBar(widgets=DEFAULT_WIDGETS, maxval=nt).start()

        for t in range(nt):
            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        data = other.read_data(t, x, y, z)
                        self.write_data(data)

            pbar.update(t)

        pbar.finish()

class PropagatorMDPSplit(QCDFormat):
    site_order = [T, X, Y, Z]
    is_gauge = False

    def __init__(self, filename):
        self.filename = filename
        self.header_format = '<60s60s60sLi10iii'
        self.endianess = '<'
        self.header_size = 60 + 60 + 60 + 14 * 4
        self.offset = None
        self.site_size = None
        self.base_size = 16 * 9 * 2

    def write_header(self, precision, nt, nx, ny, nz):
        with open(self.filename, 'wb') as f:
            self.site_size = self.base_size * (4 if precision == 'f' else 8)
            now = NOW.isoformat()
            data = struct.pack(self.header_format.encode('utf-8'),
                               'File Type: MDP FIELD'.encode('utf-8'),
                               'none'.encode('utf-8'), #self.filename.encode('utf-8'),
                               now.encode('utf-8'),
                               1325884739, 4, nt, nx, ny, nz, 0, 0, 0, 0, 0, 0,
                               self.site_size, nt * nx * ny * nz)
            f.write(data)
            self.lattice_size = (nt, nx, ny, nz)
            self.precision = precision
            self.offset = f.tell()

    def write_data(self, data, target_precision=None):
        if len(data) != self.base_size:
            raise RuntimeError("invalid data size")
        return self.file.write(self.pack(data))

    def convert_from(self, other, target_precision=None):
        precision, nt, nx, ny, nz = other.read_header()
        print(f"  (precision: {precision}, size: {nt}x{nx}x{ny}x{nz})")
        pbar = ProgressBar(widgets=DEFAULT_WIDGETS, maxval=nt).start()

        for t in range(nt):
            slice = PropagatorMDP(self.filename.replace('.split.mdp', f'.t{t:04d}.mdp'))
            slice.write_header(target_precision or precision, 1, nx, ny, nz)
            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        data = other.read_data(t, x, y, z)
                        slice.write_data(data)
            slice.close()
            pbar.update(t)
        pbar.finish()

class GaugeILDG(QCDFormat):
    def __init__(self, filename, lfn='unknown'):
        super().__init__()
        self.filename = filename
        self.endianess = '>'
        self.lfn = lfn
        self.field = 'su3gauge'
        self.base_size = 4 * 9 * 2

    def read_header(self):
        self.lime = Lime(self.filename, 'r')
        self.file = self.lime.file
        for name, stream, size in self.lime:
            if name == 'ildg-binary-data':
                self.offset = stream.tell()

        for name, stream, size in self.lime:
            if name == 'ildg-format':
                data = stream.read(size).rstrip(b'\0')
                dxml = Lime.xml_parser(data)
                field = dxml("field")
                if field != self.field:
                    raise IOError('not a lime GaugeILDG')
                precision = int(dxml("precision"))
                nt = int(dxml("lt"))
                nx = int(dxml("lx"))
                ny = int(dxml("ly"))
                nz = int(dxml("lz"))
                if precision == 32:
                    self.precision = 'f'
                    self.site_size = self.base_size * 4
                elif precision == 64:
                    self.precision = 'd'
                    self.site_size = self.base_size * 8
                else:
                    raise IOError("unable to determine input precision")
                self.lattice_size = (nt, nx, ny, nz)
                return self.precision, nt, nx, ny, nz

        raise IOError("file is not in lime format")

    def write_header(self, precision, nt, nx, ny, nz):
        self.precision = precision
        self.site_size = 4 * 2 * 9 * (4 if precision == 'f' else 8)
        self.lattice_size = (nt, nx, ny, nz)
        self.lime = Lime(self.filename, 'w')
        self.file = self.lime.file
        precision_value = 32 if precision == 'f' else 64
        d = dict(field='su3gauge', version='1.0', precision=precision_value, lx=nx, ly=ny, lz=nz, lt=nt)
        data = f"""<?xml version = "1.0" encoding = "UTF-8"?>
<ildgFormat>
<version>{d['version']}</version>
<field>{d['field']}</field>
<precision>{d['precision']}</precision>
<lx>{d['lx']}</lx><ly>{d['ly']}</ly><lz>{d['lz']}</lz><lt>{d['lt']}</lt>
</ildgFormat>"""
        self.lime.write('ildg-format', data)

    def read_data(self, t, x, y, z):
        nt, nx, ny, nz = self.lattice_size
        i = self.offset + (x + nx * (y + ny * (z + nz * t))) * self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        return self.unpack(data)

    def write_data(self, data, target_precision=None):
        if len(data) != self.base_size:
            raise RuntimeError("invalid data size")
        return self.file.write(self.pack(data))

    def convert_from(self, other, target_precision=None):
        precision, nt, nx, ny, nz = other.read_header()
        print(f"  (precision: {precision}, size: {nt}x{nx}x{ny}x{nz})")
        self.write_header(target_precision or precision, nt, nx, ny, nz)
        pbar = ProgressBar(widgets=DEFAULT_WIDGETS, maxval=nt).start()

        def reader():
            for t in range(nt):
                for z in range(nz):
                    for y in range(ny):
                        for x in range(nx):
                            data = other.read_data(t, x, y, z)
                            yield self.pack(data)
                pbar.update(t)

        self.lime.write('ildg-binary-data', reader(), nt * nx * ny * nz * self.site_size)
        self.lime.write('ildg-data-LFN', self.lfn)
        self.lime.close()
        pbar.finish()

class GaugeSCIDAC(QCDFormat):
    def __init__(self, filename):
        super().__init__(None, '>', True)
        self.filename = filename
        self.base_size = 4 * 9 * 2  # Base size for site (float precision)
        self.lattice_size = None

    def read_header(self):
        self.lime = Lime(self.filename, 'r')
        self.file = self.lime.file

        for name, stream, size in self.lime:
            if name == 'scidac-binary-data':
                self.offset = stream.tell()

        for name, stream, size in self.lime:
            if name == 'scidac-private-file-xml':
                data = stream.read(size).rstrip(b'\0')
                dxml = Lime.xml_parser(data)
                dims = dxml("dims").strip().split()
                nt, nx, ny, nz = map(int, [dims[3], dims[0], dims[1], dims[2]])
                self.lattice_size = (nt, nx, ny, nz)

        for name, stream, size in self.lime:
            if name == 'scidac-private-record-xml':
                data = stream.read(size).rstrip(b'\0')
                dxml = Lime.xml_parser(data)
                precision = dxml("precision").lower()
                if precision == 'f':
                    self.precision = 'f'
                    self.site_size = self.base_size * 4
                elif precision == 'd':
                    self.precision = 'd'
                    self.site_size = self.base_size * 8
                else:
                    raise IOError("Unable to determine input precision")

        if self.lattice_size and self.precision:
            return self.precision, *self.lattice_size
        raise IOError("File is not in lime format")

    def read_data(self, t, x, y, z):
        nt, nx, ny, nz = self.lattice_size
        i = self.offset + (x + nx * (y + ny * (z + nz * t))) * self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        return self.unpack(data)

class PropagatorSCIDAC(QCDFormat):
    is_gauge = False

    def __init__(self, filename):
        super().__init__(None, '>', False)
        self.filename = filename
        self.base_size = 16 * 9 * 2

class GaugeMILC(QCDFormat):
    def __init__(self, filename):
        super().__init__(None, '<', True)
        self.filename = filename
        self.base_size = 4 * 9 * 2  # Base size for site (float precision)
        self.header_format = '<i4i64siii'  # may change
        self.header_size = 96
        self.offset = None
        self.site_size = None

    def read_header(self):
        with open(self.filename, 'rb') as self.file:
            header = self.file.read(self.header_size)
            for self.header_format in ('<i4i64siii', '>i4i64siii'):
                self.endianess = self.header_format[0]
                items = struct.unpack(self.header_format, header)

                if items[0] == 20103:
                    nt, nx, ny, nz = items[4], items[1], items[2], items[3]
                    self.site_size = (os.path.getsize(self.filename) - 96) // (nt * nx * ny * nz)
                    self.lattice_size = (nt, nx, ny, nz)

                    if self.site_size == self.base_size * 4:
                        self.precision = 'f'
                    elif self.site_size == self.base_size * 8:
                        self.precision = 'd'
                    else:
                        raise IOError("File not in GaugeMILC format")

                    self.offset = self.file.tell()
                    return self.precision, nt, nx, ny, nz

        raise IOError("File is not in lime format")

    def read_data(self, t, x, y, z):
        nt, nx, ny, nz = self.lattice_size
        i = self.offset + (x + nx * (y + ny * (z + nz * t))) * self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        return self.unpack(data)

class GaugeNERSC(QCDFormat):
    def __init__(self, filename):
        super().__init__(None, '>', True)
        self.filename = filename
        self.offset = None
        self.site_size = None
        self.base_size = 4 * 9 * 2

    def read_header(self):
        with open(self.filename, 'rb') as file:
            header = file.read(100000)
            self.offset = header.find(b'END_HEADER\n') + 11
            if self.offset < 20:
                raise IOError('Not in NERSC format')

            lines = header[:self.offset - 1].split(b'\n')[1:-2]
            info = dict([[x.strip() for x in item.split(b' = ', 1)] for item in lines])

            nx, ny, nz, nt = map(int, [info[b'DIMENSION_1'], info[b'DIMENSION_2'], info[b'DIMENSION_3'], info[b'DIMENSION_4']])

            self.endianess = '<' if info[b'FLOATING_POINT'].endswith(b'SMALL') else '>'

            if info[b'DATATYPE'] == b'4D_SU3_GAUGE_3x3':
                self.reunitarize = False
            elif info[b'DATATYPE'] == b'4D_SU3_GAUGE':
                self.reunitarize = True
                self.base_size = 4 * 6 * 2
            else:
                raise IOError("Not in a known NERSC format")

            if info[b'FLOATING_POINT'].startswith(b'IEEE32'):
                self.precision = 'f'
                self.site_size = self.base_size * 4
            elif info[b'FLOATING_POINT'].startswith(b'IEEE64'):
                self.precision = 'd'
                self.site_size = self.base_size * 8
            else:
                raise IOError("Unable to determine input precision")

            self.lattice_size = (nt, nx, ny, nz)
            return self.precision, nt, nx, ny, nz

    def read_data(self, t, x, y, z):
        nt, nx, ny, nz = self.lattice_size
        index = self.offset + (x + nx * (y + ny * (z + nz * t))) * self.site_size
        with open(self.filename, 'rb') as file:
            file.seek(index)
            data = file.read(self.site_size)
            items = self.unpack(data)

            if self.reunitarize:
                items = [item for i in range(4) for item in reunitarize(items[i * 12:(i + 1) * 12])]

            return items

OPTIONS = {
    'mdp': (GaugeMDP, GaugeMDP, GaugeMILC, GaugeNERSC, GaugeILDG, GaugeSCIDAC),
    'ildg': (GaugeILDG, GaugeILDG, GaugeMILC, GaugeNERSC, GaugeMDP, GaugeSCIDAC),
    'prop.mdp': (PropagatorMDP, PropagatorMDP, PropagatorSCIDAC),
    'prop.ildg': (PropagatorSCIDAC, PropagatorSCIDAC, PropagatorMDP),
    'split.mdp': (GaugeMDPSplit, GaugeMDP, GaugeMILC, GaugeNERSC, GaugeILDG, GaugeSCIDAC),
    'split.prop.mdp': (PropagatorMDPSplit, PropagatorMDP, PropagatorSCIDAC),
}

def universal_converter(path, target, precision):
    filenames = [f for f in glob.glob(path) if not f.endswith(f'.{target}')]
    if not filenames:
        raise RuntimeError("No files to be converted")

    processed = set()
    messages = []
    option = OPTIONS[target]

    for filename in filenames:
        for formatter in option[1:]:
            messages.append(f'Trying to convert {filename} ({formatter})')
            try:
                ofilename = f'{filename}.{target}'
                option[0](ofilename).convert_from(formatter(filename), precision)
                processed.add(filename)
                break
            except Exception as e:
                messages.append(f'Unable to convert:\n{traceback.format_exc()}')

        if filename not in processed:
            notify('\n'.join(messages))
            sys.exit(1)

def md5_for_large_file(filename, block_size=2**20):
    """Computes MD5 hash for a large file in blocks to avoid memory overload."""
    if not os.path.exists(filename):
        return None
    with open(filename, 'rb') as f:
        md5 = hashlib.md5()
        while True:
            data = f.read(block_size)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()

def get_list(url):
    """Fetches the file list from the provided URL."""
    try:
        response = urllib.urlopen(url)
        json_data = response.read()
        files = eval(json_data)['files']
        return files
    except Exception:
        return None

def download(files, target_folder, options):
    """Downloads files from URLs specified in the 'files' list."""
    notify(f"Total files to download: {len(files)}")
    for k, f in enumerate(files):
        path = f['filename']
        basename = os.path.basename(path)
        target_name = os.path.join(target_folder, basename)

        if not os.path.exists(target_name) or os.path.getsize(target_name) != f['size']:
            input_data = None
            widgets = [basename, Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]

            while not input_data:
                notify(f"Downloading {basename}")
                try:
                    input_data = urllib.urlopen(f['link'])
                except IOError:
                    notify(f"Failure to download {f['link']}")
                    sys.exit(1)

                length = int(input_data.info().get('Content-Length', f['size']))
                if not input_data:
                    notify(f"Unable to retrieve {basename}, retrying in 5 minutes")
                    time.sleep(5 * 60)

            if not options.quiet:
                pbar = ProgressBar(widgets=widgets, maxval=length).start()

            with open(target_name, 'wb') as output:
                i = 0
                while True:
                    data = input_data.read(MAXBYTES)
                    if not data:
                        break
                    output.write(data)
                    i += len(data)
                    if not options.quiet:
                        pbar.update(i)

            input_data.close()

            if not options.quiet:
                pbar.finish()

            notify(f"Completed downloads: {k + 1}/{len(files)}")
        else:
            notify(f"Skipping file {basename} (already present)")

def ftp_download(source, target_folder, username, password):
    """
    Connects to an FTP server and downloads all files to the specified folder.
    
    :param source: FTP server address (e.g., "ftp.example.com")
    :param target_folder: Local folder to save downloaded files
    :param username: FTP username
    :param password: FTP password
    """

    # Ensure the download folder exists
    os.makedirs(target_folder, exist_ok=True)

    try:
        # Connect to FTP server
        ftp = FTP(source)
        ftp.login(username, password)

        # List files in the current directory
        files = ftp.nlst()

        notify(f"Total files to download: {len(files)}")

        for k, filename in enumerate(files):
            target_path = os.path.join(target_folder, os.path.basename(filename))

            # Get file size from FTP
            try:
                size = ftp.size(filename)
            except Exception:
                notify(f"Skipping {filename}: Unable to retrieve size.")
                continue

            # Skip if file exists with the correct size
            if os.path.exists(target_path) and os.path.getsize(target_path) == size:
                notify(f"Skipping {filename} (already downloaded)")
                continue

            notify(f"Downloading {filename}...")

            with open(target_path, 'wb') as file:
                if size:
                    pbar = ProgressBar(widgets=[filename, Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()], maxval=size).start()
                
                def write_chunk(data):
                    file.write(data)
                    if size:
                        pbar.update(file.tell())

                try:
                    ftp.retrbinary(f'RETR {filename}', write_chunk)
                except Exception as e:
                    notify(f"Failed to download {filename}: {e}")
                    os.remove(target_path)  # Remove incomplete file
                    continue

                if size:
                    pbar.finish()

            notify(f"Downloaded: {filename} ({k + 1}/{len(files)})")

        ftp.quit()
        notify("All files downloaded successfully.")

    except Exception as e:
        notify(f"Error: {e}")
        sys.exit(1)

def test_conversions():
    """Test various conversions."""
    try:
        passed = False
        test_lime()
        notify("Testing GaugeCold -> gauge_cold_1.mdp")
        GaugeMDP('gauge_cold_1.mdp').convert_from(GaugeCold(4, 4, 4, 4))
        notify("Testing gauge_cold_1.mdp -> gauge_cold_1.ildg")
        GaugeILDG('gauge_cold_1.ildg').convert_from(GaugeMDP('gauge_cold_1.mdp'))
        notify("Testing gauge_cold_1.ildg -> gauge_cold_2.mdp")
        GaugeMDP('gauge_cold_2.mdp').convert_from(GaugeILDG('gauge_cold_1.ildg'))
        notify("Comparing gauge_cold_1.mdp with gauge_cold_2.mdp")
        assert open('gauge_cold_1.mdp', 'rb').read() == open('gauge_cold_2.mdp', 'rb').read()
        GaugeMDP('gauge_cold_3.mdp').convert_from(GaugeMDP('gauge_cold_2.mdp'))
        notify("Comparing gauge_cold_1.mdp with gauge_cold_3.mdp")
        assert open('gauge_cold_1.mdp', 'rb').read() == open('gauge_cold_3.mdp', 'rb').read()
        GaugeILDG('gauge_cold_2.ildg').convert_from(GaugeILDG('gauge_cold_1.ildg'))
        notify("Comparing gauge_cold_1.ildg with gauge_cold_2.ildg")
        assert open('gauge_cold_1.ildg', 'rb').read() == open('gauge_cold_2.ildg', 'rb').read()
        GaugeMDP('gauge_cold_4.mdp').convert_from(GaugeNERSC('samples/demo.nersc'))
        GaugeILDG('gauge_cold_4.ildg').convert_from(GaugeNERSC('samples/demo.nersc'))
        GaugeMDP('gauge_cold_5.mdp').convert_from(GaugeILDG('gauge_cold_4.ildg'))
        notify("Comparing gauge_cold_4.mdp with gauge_cold_5.mdp")
        assert open('gauge_cold_4.mdp','rb').read() == open('gauge_cold_5.mdp','rb').read()
        passed = True
    except Exception as e:
        notify(f"Tests failed:\n{traceback.format_exc()}")
    finally:
        notify("Cleaning up tmp files")
        os.system('rm gauge_cold_?.*')

    if not passed:
        sys.exit(1)

class DummyProgressBar(object):
    """Used to avoid the ProgressBar on -n, --noprogressbar."""
    def __init__(self, *a, **b):
        self.maxval = b['maxval']

    def update(self, i):
        notify(f"...{i}/{self.maxval}")

    def start(self):
        return self

    def finish(self):
        notify("...done!")

def main():
    """Main program to handle file download, conversion, and tests."""
    print(LOGO)
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-q", "--quiet", dest='quiet', action='store_true', default=False, help='No progress bars')
    parser.add_option("-d", "--destination", dest='destination', default=None, help="Destination folder")
    parser.add_option("-c", "--convert", dest='convert', default=False, help="Convert a field to format (mdp, ildg)")
    parser.add_option("-4", "--float", dest='float_precision', default=False, action='store_true', help="Convert to float precision")
    parser.add_option("-8", "--double", dest='double_precision', default=False, action='store_true', help="Convert to double precision")
    parser.add_option("-t", "--tests", dest='tests', default=False, action='store_true', help="Run some tests")
    parser.add_option("-n", "--noprogressbar", dest='noprogressbar', default=False, action='store_true', help="Disable progress bar")
    (options, args) = parser.parse_args()

    ### disable progress bar if necessary
    if options.noprogressbar:
        global ProgressBar
        ProgressBar = DummyProgressBar

    ### run tests if asked
    if options.tests:
        test_conversions()
        sys.exit(0)

    ### if no argument source passed, print usage
    try:
        options.source = args[0]
    except IndexError:
        print(USAGE)
        sys.exit(1)

    ### download data (http, https, ftp, sftp) or not
    if options.source.startswith(('http://', 'https://')):
        files = get_list(options.source)
        if files is None:
            notify("Unable to connect")
            sys.exit(0)
        else:
            regex = re.compile('pattern = (?P<pattern>[^&]*)')
            pattern = regex.search(options.source).group('pattern')
            target_folder = options.destination or urllib.unquote(pattern)
            notify(f"Target folder: {target_folder}")
            if not os.path.exists(target_folder):
                os.mkdir(target_folder)
            download(files, target_folder, options)
        conversion_path = os.path.join(target_folder, pattern.replace('nnnnn', '*'))
    elif options.source.startswith(('ftp://', 'sftp://')):
        username = input('Username: ')
        password = input('Password: ')
        if not os.path.exists(target_folder):
            os.mkdir(target_folder)
        ftp_download(options.source, target_folder, username, password)
    else:
        conversion_path = options.source

    ### if options.source = 'gauge.cold.TxXxYxZ' make it
    if options.source.startswith('gauge.cold') and not os.path.exists(options.source):
        size = [int(x) for x in options.source[11:].split('x')]
        GaugeMDP(options.source).convert_from(GaugeCold(*size))

    ### if conversion required use the universal converter
    if options.convert:
        notify(f"Converting: {conversion_path} -> {conversion_path}.{options.convert}")
        precision = 'f' if options.float_precision else 'd' if options.double_precision else None
        universal_converter(conversion_path, options.convert, precision)

if __name__ == '__main__':
    main()
