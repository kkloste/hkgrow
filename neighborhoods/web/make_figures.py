
project = 'Neighborhood Communities'
page = 'Online Figures'
title = '%s | %s'%(project,page)
authors = 'David F. Gleich and C. Seshadhri'

HEADER="""
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
        "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en">
<head>
        
       <title>%(title)s</title>

       <link rel="stylesheet" href="css/lightbox.css" type="text/css" media="screen" />
        
        <script src="js/prototype.js" type="text/javascript"></script>
        <script src="js/scriptaculous.js?load=effects,builder" type="text/javascript"></script>
        <script src="js/lightbox.js" type="text/javascript"></script>
</head>
<body>
"""%(locals())

FOOTER="""
<emph>Back to <a href=".">%s project home</a></emph>
</body>
</html>
"""%(project)

print HEADER

print """
<h1>%(project)s</h1>
<h2>%(page)s</h2>
<emph>%(authors)s</emph>
"""%(locals())

def print_set(figpre,title):
    graphs = [['Collaboration','ca-AstroPh-cc','email-Enron-cc',
                'cond-mat-2005-fix-cc','arxiv-ubc','dblp-cc',
                'hollywood-2009-cc'],
              ['Social','Penn94','anony-interactions-oneyearA-cc',
                'networkA-anonymized','soc-LiveJournal1'],
              ['Technological','oregon2_010526','p2p-Gnutella25',
                'as-22july06', 'itdk0304-cc'],
              ['Web','web-Google'],
              ['Models','rand-ff-25000-0.4','rand-ff-25000-0.49']]
    print "<h3>%s</h3>"%(title)
    print '<a name="%s" />'%(figpre)
    for types in graphs:
        curtype = types[0]
        print "<p>%s</p>"%(curtype)
        print "<ul>"
        for g in types[1:]:
            print '<a href="%s-%s.png" rel="lightbox[%s]" title="%s">%s</a>'%(
                figpre, g, figpre, g, g)
        print "</ul>"
    

print_set('neigh','Neigborhood communities');
print_set('nvsncp','Neighborhoods compared to Personalized PageRank');
print_set('core','k-core communities');
print_set('pprgrow','PPR communities from locally minimal seeds');


print FOOTER

