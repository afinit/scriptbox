/^[ATCG]/ { l+=length }
/^>/ {  if( l > 0 ){
          if( l >= 1000 ){
            c4+=1;
            t4+=l;
          }
          if( l < 1000 ){
            c+=1;
            t+=l;
          }
          c2+=1;
          t2+=l;
          l=0;
        }
     }

END { if( l < 1000 ){
        c+=1;
        t+=l;
      }
      if( l >= 1000 ){
        c4+=1;
        t4+=l;
      }
      c2+=1;
      t2+=l;
      
      if( c4 > 0 ){
        print "tot contigs >= 1000:\t" c4;
        print "avg length >= 1000 :\t" t4/c4;
        print "tot bp's >= 1000   :\t" t4;
        print " ";
      }
      if( c > 0 ){
        print "tot contigs under 1000:\t" c;
        print "avg length under 1000 :\t" t/c;
        print "tot bp's under 1000   :\t" t;
        print " ";
      }
      print "tot contigs:\t" c2;
      print "avg length:\t" t2/c2;
      print "tot bp's:\t" t2; 
    }
