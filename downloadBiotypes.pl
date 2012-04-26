use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all("/data/ensembl/registry/ensembl_registry.conf");

my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

my $dba = $registry->get_DBAdaptor("Human", "core");

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

my $simple_feature_adaptor = $dba->get_SimpleFeatureAdaptor();

@slices = @{ $slice_adaptor->fetch_all('chromosome') }; 

#$rca = $dba->get_RepeatConsensusAdaptor;

## For getting the type of repeats in the genome uncomment the following code <repeatsTypes>
#my ($sth, @types, $repeat_type);
#$dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
      #-user   => 'anonymous',
      #-dbname => 'homo_sapiens_core_64_37',
      #-host   => 'ensembldb.ensembl.org',
      #-driver => 'mysql',
      #-port => '5306'
      #);

#$sth = $dbc->prepare("SELECT DISTINCT repeat_type FROM repeat_consensus");
#$sth->execute();
#$sth->bind_columns(\$repeat_type);
#while ($sth->fetch()) { push @types, $repeat_type }

#foreach my $type (@types){
	#print "$type\n";
#}

## </repeatsTypes>

#foreach my $slice (@slices){
	#my $genes = $slice->get_all_Genes();

	#foreach my $gene (@$genes) {
		 #printf( "%s\t%s\n", $gene->stable_id(),$gene->biotype());
	#}
#}


my $dba_mouse = $registry->get_DBAdaptor("Mouse", "core");

my $slice_adaptor_mouse = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );

my $simple_feature_adaptor_mouse = $dba_mouse->get_SimpleFeatureAdaptor();

@slices_mouse = @{ $slice_adaptor_mouse->fetch_all('chromosome') }; 

foreach my $slice (@slices){
	my $genes = $slice->get_all_Genes();

	foreach my $gene (@$genes) {
		 printf( "%s\t%s\n", $gene->stable_id(),$gene->biotype());
	}
}
