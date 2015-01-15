BEGIN { m=10000 }
/^@$/ {getline; c+=1; len=length($0); s+=len; getline; getline;
        if( m>len ){
          m=len;
        }
      }
END { print "read count:\t\t" c
      print "avg length:\t" s/c
      print "min length:\t" m }
