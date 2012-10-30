#!/usr/ensembl/bin/perl -w
#
# Written by Jan-Hinnerk Vogel
#
# Copyright GRL/EBI 2004
# You may distribute this script under the same terms as perl itself
#
#
#



use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;


my %opt = (
            dbhost => 'ensembldb.ensembl.org',
            dbname => 'mysql',
            dbuser => 'anonymous',
            dbpass => undef,
            dbport => '5306',
          );

&GetOptions(
            'dbname=s' => \$opt{dbname},
            'dbuser=s' => \$opt{dbuser},
            'dbhost=s' => \$opt{dbhost},
            'dbport=s' => \$opt{dbport},
            'dbpass=s' => \$opt{dbpass},
            'verbose+' => \$opt{verbose},
           );


my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
                                                 '-dbname' => $opt{dbname},
                                                 '-host' => $opt{dbhost},
                                                 '-user' => $opt{dbuser},
                                                 '-port' => $opt{dbport},
                                                 '-pass' => $opt{dbpass}
                                                );




my $query =qq {
	show databases like '%core%'
              };


my $sth = $dbc->prepare($query);

$sth->execute;

my @aref = @{$sth->fetchall_arrayref()};

for(@aref){
  my @row = @$_;
  print join("\t", @row) ."\n" ;
}


