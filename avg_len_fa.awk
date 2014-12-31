/^>/ { sub(/^.*cov_/, "" ); sub(/_.*$/,""); t+=$0; c+=1  }
/^[ATCG]/ { l += length }
END { print "total cont:\t" c
      print "tot bp's  :\t" l
      print "avg length:\t" l/c
      print "avg cov:\t" t/c }
