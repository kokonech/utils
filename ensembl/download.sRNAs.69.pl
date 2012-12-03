use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all("/data/ensembl/registry/ensembl_registry.69.conf");

my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

my $dba = $registry->get_DBAdaptor("Human", "core");

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

my $simple_feature_adaptor = $dba->get_SimpleFeatureAdaptor();

@slices = @{ $slice_adaptor->fetch_all('chromosome') }; 

my @sRNAs_biotypes = qw(miRNA misc_RNA  Mt_tRNA rRNA snoRNA  snRNA); # taken from http://www.ensembl.org/info/docs/Doxygen/core-api/BiotypeMapper_8pm_source.html

foreach my $slice (@slices){
	my $genes = $slice->get_all_Genes();

	foreach my $gene (@$genes) {
	    
        if ($gene->biotype() ~~ @sRNAs_biotypes){
            printf( "%s\t%s\t%s\t%d\n", $gene->stable_id(),$gene->external_name(),$gene->biotype(),$gene->length());
        }

	}
}
