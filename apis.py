import sys
import json
import time
import re

# Python 2/3 adaptability

from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError

import requests


class EnsemblRestClient(object):
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = Request(self.server + endpoint, headers=hdrs)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data
    def get_overlap(self, region):
        chromosome, loc = [None, None]
        try:
            chromosome, loc = region.split(":")
        except:
            chromosome, loc, blah = region.split(":")
        reg = "{0}:{1}-{1}".format(chromosome, loc)
        overlap = self.perform_rest_action(
            endpoint='/overlap/region/human/{}'.format(reg),
            params={"feature":"variation"}
        )
        return overlap
    def get_effects(self, species, symbol):
        effects = self.perform_rest_action(
            endpoint='/vep/{0}/id/{1}'.format(species, symbol)
            #params={'object_type': 'gene'}
        )
        
        """if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants"""
        return effects


class GwasRestClient():
    def __init__(self, server="https://www.ebi.ac.uk/gwas/rest/api", reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
    
    def check_status(self, jobId):
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        ids = None
        ids = self.perform_rest_action(endpoint="/idmapping/status/" + jobId)
        if (not "results" in ids.keys()):
            ids = self.check_status(jobId)
        return ids
    
    
    def perform_rest_action(self, endpoint, hdrs=None, form_data=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'
        
        
        if form_data != None:
            form_data = urlencode(form_data)
            form_data = bytes(form_data, 'ascii')
        
        data = None
        
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        try:
            request = Request(self.server + endpoint, data=form_data)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_assoc(self, efoTrait):
        assoc = self.perform_rest_action(endpoint="/efoTraits/{0}/associations".format(efoTrait))
        if assoc:
            return assoc
        else:
            print("failed to get associations")
    
    #transcript_ids are comma-separated
    def id_map(self, transcript_ids):
        data = {
            "ids": transcript_ids,
            "from": "Ensembl_Transcript",
            "to": "UniProtKB"
        }
        
        jobId = self.perform_rest_action(endpoint="/idmapping/run", form_data=data)
        if jobId:
            jobId = jobId["jobId"]
            ids = self.check_status(jobId)
            if ids:
                return ids
        else:
            print("ID mapping failed")
    

class UniProtRestClient():
    def __init__(self, server="https://rest.uniprot.org", reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
    
    def get_all_pages(self, url):
        response = requests.get(url)
        hdrs = response.headers
        next_page_data = self.get_all_pages(hdrs["Link"]) if "Link" in hdrs.keys() else ""
        data = response.text
        return data + next_page_data


    def check_status(self, jobId):
        ids = self.perform_rest_action(endpoint="/idmapping/status/" + jobId)
        while (not "results" in ids.keys()):
            if self.req_count >= self.reqs_per_sec:
                delta = time.time() - self.last_req
                if delta < 1:
                    time.sleep(1 - delta)
                self.last_req = time.time()
                self.req_count = 0
            ids = self.perform_rest_action(endpoint="/idmapping/status/" + jobId)
        #print(jobId)
        #"https://rest.uniprot.org/idmapping/uniprotkb/results/8cf41d545f303e03d20460d27b59dcd1d674bcec?format=json&size=500"
        
        ids_url = f"{self.server}/idmapping/uniprotkb/results/stream/{jobId}?format=json&size=500"
        ids = requests.get(ids_url).json()
        url = "{0}/idmapping/uniprotkb/results/{1}?compressed=false&format=fasta&query=%28reviewed%3Atrue%29&size=500".format(self.server, jobId)
        all_fastas = self.get_all_pages(url)
        fasta_list = re.split(r'\n(?=>)', all_fastas)
        [fasta for fasta in fasta_list]
        return fasta_list, ids
    
    
    def perform_rest_action(self, endpoint, form_data=None):
        if form_data != None:
            form_data = urlencode(form_data)
            form_data = bytes(form_data, 'ascii')
        
        data = None
        
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        try:
            request = Request(self.server + endpoint, data=form_data)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    
    #transcript_ids are comma-separated
    def id_map(self, transcript_ids):
        data = {
            "ids": transcript_ids,
            "from": "Ensembl_Transcript",
            "to": "UniProtKB"
        }
        
        jobId = self.perform_rest_action(endpoint="/idmapping/run", form_data=data)
        if jobId:
            jobId = jobId["jobId"]
        else:
            print("ID mapping failed")
            return
        fasta_list, ids = self.check_status(jobId)
        if fasta_list:
            return fasta_list, ids
     