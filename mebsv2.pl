#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin '$Bin';

# General script to score biogeochem cycles in both genomic and metagenomic data.
# B Contreras-Moreira, V de Anda 2018
# bcontreras@eead.csic.es , valdeanda@ciencias.unam.mx 

my $HMMSEARCHEXE = 'hmmsearch'; # please edit if not in path
my $CONFIGDIR   = $Bin.'/config/';
my $CONFIGFILE  = $CONFIGDIR . 'config.txt';
my $VALIDEXT    = '.faa';
my $VALIDENT    = '.tab';
my $HMMOUTEXT   = '.hmmsearch.tab'; # default extension for hmmsearch tbl output
my @validMSL    = qw(30 60 100 150 200 250 300);
my ($INP_help,$INP_folder,$INP_cycles,$INP_type, $INP_db, $INP_comp) = (0,'',0,'',0,'');
my $db         = 'Pfam-A.hmm';
GetOptions
(
  'help|h|?'    => \$INP_help,
  'input|in=s'  => \$INP_folder,
  'type|t=s'    => \$INP_type,
  'comp|mc'     => \$INP_comp,
  'db|d'        => \$INP_db,
);


if (-t STDIN && ($INP_help || $INP_folder eq '' || $INP_type eq '' || $INP_db eq '', ) && !$INP_cycles)
{
  die<<EODOC;

  Program to compute MEBS for a set of genomic/metagenomic FASTA files in input folder.
  
  usage: $0 [options] 

   -help    Brief help message
   
   -input   Folder containing FASTA peptide files ($VALIDEXT)             (required)

   -type    Nature of input sequences, either 'genomic' or 'metagenomic'  (required)
   
   -db      Database to scann your input files                (required, default $db)
   

   -comp    Compute the metabolic completeness                            (optional)


EODOC
}

## 1) Checking binaries
if(!$HMMSEARCHEXE )
{
  die "#ERROR:  hmmsearch not found, please install or set \$HMMSEARCHEXE correctly\n";
}

if(!$INP_db)

{
  die "ERROR: db not found, please install manually see \n
      ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz";
}


## 2) Checking parameters
my (@valid_infiles, @cycles, @config, @paths, @completeness);
my (@MSL,  %pathways);
my ($c,$f,$path,$cycle,$msl,$score,$comp,$pw);
my ($hmmsearchfile,$entropyfile,$scorefile,$infile,$pfam2keggfile);


# Read config file
open(CONFIG,$CONFIGFILE) || die "# ERROR: cannot read $CONFIGFILE\n";
while(my $line = <CONFIG>)
{
  #Cycle   Path    Comple  Input Genes     Input Genomes   Domains AUC     Score(FD..
  #sulfur  cycles/sulfur/  cycles/sulfur/pfam2kegg.tab     152     161     112    ..

  @config = split(/\t/,$line);
  if($config[0] =~ /^Cycle/)
  {
  }
  else
  {
    push(@cycles, $config[0]);
    push(@paths, $config[1]);
    push(@completeness, $config[2]);  
   }
}  

close(CONFIG);


print   warn "# $0 -input $INP_folder -type $INP_type  -comp $INP_comp\n\n";
 
# check required sequence type
if(!$INP_type || ($INP_type ne 'genomic' && $INP_type ne 'metagenomic'))
{
  die "# ERROR : type of input must be indicated; valid options are [genomic|metagenomic]\n";
}

# check required sequence folder
if(!$INP_folder)
{
  die "# ERROR : need valid -input folder\n";
}
else
{
  opendir(DIR,$INP_folder) || die "# ERROR: cannot list $INP_folder\n";
  @valid_infiles = grep{/$VALIDEXT$/} readdir(DIR);
  closedir(DIR);
  if(scalar(@valid_infiles) == 0)
  {
    die "# ERROR: cannot find files with extension $VALIDEXT in folder $INP_folder\n";
  }
  elsif($INP_type eq 'metagenomic')
  {
    # compute Mean Size Length for this metagenomic sequence set
    warn "# Computing Mean Size Length (MSL) ...\n";

    my ($c,$nseq,$mean,$len,$cutoff,@MSLcutoffs);
    for(my $bin=0;$bin<scalar(@validMSL)-1;$bin++)
    {
      $cutoff = $validMSL[$bin] + (($validMSL[$bin+1]-$validMSL[$bin])/2);
      push(@MSLcutoffs,$cutoff);#print "$validMSL[$bin] $cutoff\n";
    }
    push(@MSLcutoffs,$validMSL[$#validMSL]); # add last MSL

    foreach my $infile (@valid_infiles)
    {
      ($nseq,$mean,$len) = (0,0,0);
      open(FAAFILE,"<","$INP_folder/$infile") || 
        die "# ERROR: cannot read files $INP_folder/$infile\n";
      while(<FAAFILE>)
      {
        if(/^>/)
        {
          $nseq++;
          $mean += $len;
          $len=0;
        }
        else
        {
          chomp;
          $len += length($_);
        }
      }
      close(FAAFILE);

      $mean = sprintf("%1.0f",$mean/$nseq);

      # find out which pre-defined MSL bin matches the estimated MSL for this sequence set
      foreach $c (0 .. $#MSLcutoffs)
      {
        $cutoff = $MSLcutoffs[$c];
        if($mean <= $cutoff)
        {
          push(@MSL,$validMSL[$c]);
          warn "# $infile MSL=$mean MSLbin=$validMSL[$c]\n";

          last;
        }
      }
    }
  } print "\n";
}


## 3) scan input sequences with selected database 

# print header
my $pathways_header = '';
foreach $c (0 .. $#cycles)
{
  print "\t$cycles[$c]";

  # print completeness header if required
  $comp = $completeness[$c];
  if($INP_comp && $comp ne "" && -s $comp)
  {
    open(COMPFILE,"<",$comp) || warn "# ERROR: cannot read $comp\n";
    while(<COMPFILE>)
    {
      #PFAM  KO  PATHWAY   PATHWAY NAME 
      #PF00890 K00394  1 Sulfite oxidation 
      if(/^PF\d+\t.*?\t(\d+)\t/)
      {
        $pathways{$cycles[$c]}{$1} = 1; 
      }
    }
    close(COMPFILE);
     
     $pathways_header ="\t$cycles[$c]_<comp>";
    foreach $pw (sort {$a<=>$b} keys(%{$pathways{$cycles[$c]}}))
    {
      $pathways_header .= "\t$cycles[$c]_$pw";
    }
  }
} print "$pathways_header\n"; 

foreach $f (0 .. $#valid_infiles)
{
  $infile = $valid_infiles[$f];

  print "$infile"; # rowname

  # compute & print scores per cycle
  foreach $c (0 .. $#cycles)
  {
    $path = $paths[$c];
    print $path; 
    $cycle = $cycles[$c];
    $comp = $completeness[$c];
    $score = 'NA';

    $hmmsearchfile = $INP_folder . '/' . $infile . '.' . $cycle . $HMMOUTEXT;
    $scorefile = $INP_folder . '/' . $infile . '.' . $cycle . '.score';
    $entropyfile = $path . 'entropies' . $VALIDENT;

    system("$HMMSEARCHEXE --cut_ga -o /dev/null --tblout $hmmsearchfile  $INP_db  $INP_folder/$infile");

    if(-s $hmmsearchfile)
    {
      if($INP_type eq 'metagenomic')
      {
        if($INP_comp && $comp ne "" && -s $comp)
        { 
          system("$Bin/scripts/pfam_score.pl -input $hmmsearchfile -entropyfile $entropyfile -size $MSL[$f] -keggmap $comp > $scorefile");
        }
        else
        {
          system("$Bin/scripts/pfam_score.pl -input $hmmsearchfile -entropyfile $entropyfile -size $MSL[$f] > $scorefile");
        }
      }
      else
      {
        if ($INP_comp && $comp ne "" && -s $comp)
        {
          system("$Bin/scripts/pfam_score.pl -input $hmmsearchfile -entropyfile $entropyfile -size real -keggmap $comp > $scorefile");
        }
        else
        {
          system("$Bin/scripts/pfam_score.pl -input $hmmsearchfile -entropyfile $entropyfile -size real > $scorefile");
        }  
      }
      
      if(-s $scorefile)
      {
        $score = -1;
        open(SCORE,"<",$scorefile) || warn "# ERROR: cannot read $scorefile\n";
        while(<SCORE>)
        {
          if(/Pfam entropy score: (\S+)/){ $score = sprintf("%1.3f",$1) } 
        }
        close(SCORE);  
      }
      else { warn "# ERROR: failed to generate $scorefile\n"  }
    }
    else { warn "# ERROR: failed to generate $hmmsearchfile\n" }

    print "\t$score";
  }

  # print completeness summary per cycle
  if ($INP_comp)
  {
    foreach $c (0 .. $#cycles)
    {
      $cycle = $cycles[$c];

      # parse score file and check whether completeness report is there
      my ($compOK,%comp_scores) = (0);
      $scorefile = $INP_folder . '/' . $infile . '.' . $cycle . '.score';
      open(COMPL,"<",$scorefile);
      while(<COMPL>)
      {   
        # path_number path_name total_domains matched %completeness matched_Pfam_domains
        #1 Sulfite oxidation   9 3 33.3  PF00890,PF12838,PF01087
        # mean pathway completeness: 37.9
        if(/^# path_number/){ $compOK = 1 }
        elsif($compOK && /^(\d+)\t.*?\t\d+\t\d+\t(\S+)/)
        {
          $comp_scores{$1} = $2;
        }
        elsif(/^# mean pathway completeness: (\S+)/)
        {
          print "\t$1"; # print mean 
        }
      }
      close(COMPL);
    
      # print completeness for all sorted pathways
      foreach $pw (sort {$a<=>$b} keys(%{$pathways{$cycle}}))
      {
        print "\t$comp_scores{$pw}";
      }  
    }
  }    

 print "\n";
}























