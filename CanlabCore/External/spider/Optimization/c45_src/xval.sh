#csh

#---------------------------------------------------------------------
# N-way cross-validation script
#---------------------------------------------------------------------
#
# invocation:
#   csh xval.sh filestem N [options for c4.5 and c4.5rules] [suffix]
#
# individual results from each block are left in
#     filestem.[rt]o*[suffix],
# averages over all blocks in
#     filestem.[rt]res[suffix]
#---------------------------------------------------------------------

#	sort the options into result suffix and control options for the programs
#	Note: for options with values, there must be no space between the option
#	name and value; e.g. "-v1", not "-v 1"

set treeopts =
set ruleopts =
set suffix =

foreach i ( $argv[3-] )
  switch ( $i )
  case "+*":
    set suffix = $i
    breaksw
  case "-v*":
  case "-c*":
    set treeopts = ($treeopts $i)
    set ruleopts = ($ruleopts $i)
    breaksw
  case "-p":
  case "-t*":
  case "-w*":
  case "-i*":
  case "-g":
  case "-s":
  case "-m*":
    set treeopts = ($treeopts $i)
    breaksw
  case "-r*":
  case "-F*":
  case "-a":
    set ruleopts = ($ruleopts $i)
    breaksw
  default:
    echo "unrecognised or inappropriate option" $i
    exit
  endsw
end

#	prepare the data for cross-validation

cat $1.data $1.test | xval-prep $2 >XDF.data
cp /dev/null XDF.test
ln $1.names XDF.names
rm $1.[rt]o[0-9]*$suffix
set junk = `wc XDF.data`
set examples = $junk[1]
set large = `expr $examples % $2`
set segsize = `expr \( $examples / $2 \) + 1`

#	perform the cross-validation trials

set i = 0
while ( $i < $2 )
  if ( $i == $large ) set segsize = `expr $examples / $2`
  cat XDF.test XDF.data | split -`expr $examples - $segsize`
  mv xaa XDF.data
  mv xab XDF.test

  c4.5 -f XDF -u $treeopts >$1.to$i$suffix
  c4.5rules -f XDF -u $ruleopts >$1.ro$i$suffix

  @ i++
end

#	remove the temporary files and summarize results

rm -f XDF.*
cat $1.to[0-9]*$suffix | grep "<<" | average >$1.tres$suffix
cat $1.ro[0-9]*$suffix | grep "<<" | average >$1.rres$suffix
