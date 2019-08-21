Roles represent:
1 step name
2 job-id
3 job create time
4 app/applet id
5 app version of applet built time
6 upstream version

Caveat: job finding use name matching and report only the last job. Unless every analysis are done on DNAnexus with its native job name, there is no guarantee that we will catch the correct execution. E.g. fNotCel1 purge_dups was done outside DNAnexus, so it would report purge_haplotig as execution even though latest results rely on purge_dups.