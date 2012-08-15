#!/usr/bin/python
'''
Created on Jul 26, 2012

@author: schudoma
'''

import os
import re
import sys
import time
import urllib2


pattern_raw = re.compile('<a href=\'[^>]+\.txt\' target=\'_new\'>View Data in Raw Text Format')
pattern_fn = re.compile('href=\'.+\.txt\'') 

base_url = 'http://bbc.botany.utoronto.ca/affydb/cgi-bin/affy_db_exprss_browser_out.cgi?'
base_params = [['datasets', 'atgenexp_plus'],
               ['output', 'ratio_avg'],
               ['agi_input', ''],
               ['pub', '']]


test_url = 'http://bbc.botany.utoronto.ca/affydb/cgi-bin/affy_db_exprss_browser_out.cgi?datasets=atgenexp_plus&output=ratio_avg&agi_input=At1g01110%0AAt1g01120&pub='

string = """
<a href='http://bbc.botany.utoronto.ca:80/ntools/temp_general/exprss_text_unsorted_30641.txt' target='_new'>View Data in Raw Text Format</a>
""".strip()

def process_url(url):
    url_obj = urllib2.urlopen(url)
    try:
        fn = get_filename(url_obj.read())
    except:
        print 'Getting filename failed for', url
        return None
    finally:
        url_obj.close()
    try:
        data = get_expression_data(fn)
    except:
        print 'Getting expression data failed for', url
        return None
    return data

def get_filename(text):
    reobj_raw = pattern_raw.search(text)
    try:
        text_raw = reobj_raw.group()
    except:
        print 'Failed raw extraction:', text
        return None
    reobj_fn = pattern_fn.search(text_raw)
    try:
        text_fn = reobj_fn.group()
    except:
        print 'Failed fn extraction:', text, text_raw
        return None
    return text_fn.lstrip('href=').strip("'")

def get_expression_data(fn):
    url = urllib2.urlopen(fn)
    data = url.read()
    url.close()
    return data

def format_agi(agi):
    return 'At%cg%s' % (agi[2], agi[4:])
    
def main(argv):    
    
    path = '/home/schudoma/ngstemp/cschudoma/Ath_Ped-0-Col-0-graft_RNA-Seq_Scheible/2012-07-10'
    fn_in = 'tair_ids.txt'
    fn_out = 'expression_data.csv'
    
    
    atg_lines = [format_agi(line.strip()) 
                 for line in open(os.path.join(path, fn_in)).readlines()]
    p = 0
    stepsize = 125
    step = stepsize
    
    try:
        # last = open('LAST_AGI.dat', 'r').read().strip()
        last = open(fn_out).readlines()[-1].split()[0]
    except:
        last = None
    if not last is None:
        p = atg_lines.index(last) + 1
    print 'P=', p
    # return None
    
    fo = open(os.path.join(path, fn_out), 'a')
    
    while p < len(atg_lines):
        print 'Current batch...', p, p+step
        batch = atg_lines[p:p+step]
        
        params = base_params
        # param format: 'At1g01110%0AAt1g01120'       
        params[2][1] = '%0A'.join(batch)
        params = ['%s=%s' % tuple(param) for param in params]        
        url = base_url + '&'.join(params)
        
        try:
            data = process_url(url)
        except:
            data = None            
            print 'Failed to process current batch.'
            # Adjusting step size ->', step / 2
            if step == 1:
                print 'Skipping id:', batch[0]
                p += 1
                step = stepsize
            else:
                print 'Adjusting step size:', step, '->', step/2
                step /= 2
            pass
        
        if not data is None:
            fo.write(data)
            # current step is still of modified stepsize (otherwise we skip over the next stepsize elements)
            p += step
            # now restore original stepsize
            step = stepsize
        
        print 'Sleeping 30s in order to confuse the server.. >;D'
        time.sleep(30)
        
        pass
    
    fo.close()
    # print count == len(atg_lines) 
        
     
        
        
    pass

def main_test(argv):
    
    
    url = urllib2.urlopen(test_url)
    
    
    text = url.read()
    url.close()
    
    print text
    print '****\n'*3,
    
    filename = get_filename(text)
    print 'URL:', filename
    data = get_expression_data(filename)
    print 'DATA:', data
    pass

if __name__ == '__main__': main(sys.argv[1:])
