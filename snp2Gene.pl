#!/usr/bin/perl

###################################################################################################################################
#
# Copyright 2017 IRD
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD
# Version 1 written by Francois Sabot
#
###################################################################################################################################

use strict;

use Data::Dumper;
use Getopt::Long;
use List::MoreUtils;

my ($snpFile,$gff, $vcf, $knet, $distance, $converter,$outFile,$help);

my $courriel="francois.sabot-at-ird.fr";
my ($nomprog) = $0 =~/([^\/]+)$/;
my $MessAbruti ="\nUsage:
\t$nomprog -i SNPfile -g gffFile -o outfile [-d distance -v vcfFile -k knetMinerFile -c IDconverterfile]

From a SNP file extracted from the GWAS analysis, will recover from the gffFile the closest gene.

The -i option is the SNPfile (mandatory)

The -g option is the corresponding GFF file (mandatory)

The -o option will specify the output file name (mandatory)

The -d option will specify the distance to 'scan', in base (standard at 1000)

The -v option will use a VCF file with the EFF/snpEff field (optional)

The -k option will allow to provide a knetMiner output (optional). If working on rice, knetminer needs also the converter from RAP-DB to MSU is GFF is not the good one. Will do it automatically, but can be overidded using the -c option

The output will be a text tabulated file such as: 
CHR     Position        SNPid   pValue  GeneID  Distance        GeneNotes	REF	ALT	QUAL	EFF	knetScore QTL	TO	Phenotype	BioProc	PMID

A negative distance in the output means in 5' of the gene, a positive one in 3' and a 0 value means within the gene.

        contact: $courriel\n\n";


#Standard values

unless (@ARGV) 
        {
        print "\nType --help for more informations\n\n";
        exit;
        }

$vcf = "null";
$knet = "null";
$converter = "/data2/projects/irigin/association/filesForGeneAnalysis/RAP-MSU_2017-08-04.txt";

GetOptions("prout|help|?|h" => \$help,   
            "i|in=s"=>\$snpFile,
            "o|out=s"=>\$outFile,
			"g|gff=s"=>\$gff,
			"v|vcf=s"=>\$vcf,
			"k|knet=s"=>\$knet,
			"c|converter=s"=>\$converter,
			"d|distance=s"=>\$distance);

if ($help)
	{
	print $MessAbruti,"\n";
	exit;
	}

#Opening files, checking right, toussa
open (my $fhIn, "<", $snpFile) or die ("\nCannot open the SNP file $snpFile:\n\n\t$!\n\nAborting...\n");
open (my $fhOut, ">", $outFile) or die("\nCannot create the output file $outFile:\n\n\t$!\n\nAborting...\n");
print $fhOut "#CHROM    POS     SNPid   pValue  GeneID  Distance GeneNotes";
print $fhOut "	REF	ALT	QUAL	ANN";
print $fhOut "	knetScore	QTL	TO	Phenotype	BioProc	PMID";

print $fhOut "\n";
#testing bedtools usability
my $bedtoolsTest = `windowBed -h 2>&1`;
chomp $bedtoolsTest;
if ($bedtoolsTest =~ m/bash:/) #No bedtools installed
{
        print "\nThis tool required BEDTools to be installed and accessible in the path to work.\nPlease contact your administrator to install BEDTools at http://bedtools.readthedocs.io/\n\n...";
        exit;
}

#READING the provided VCF file

my %vcfInfos;
if ($vcf ne "null")
{
	&vcfReader;
}

#Reading knetMiner infos
my %knetInfos;
if ($knet ne "null")
{
	&knetReader;
}

#print Dumper(\%knetInfos);
#exit;

#preparing the file as a BED file
my $tmpFile = "/tmp/tmpSNP2Gene_".(int(rand(100000))).".BED";
open (my $tmpFh, ">", $tmpFile) or die ("\nCannot create tmp files:$!\n");

while (my $line = <$fhIn>)
{
        chomp $line;
        next if $line =~ m/^SNP/;
        next if $line =~ m/^$/;
        my @data = split /\s+/, $line;
        my $outline = $data[1]."\t".$data[2]."\t".$data[2]."\t".$data[0]."\t".$data[3]."\n";
        
        print $tmpFh $outline;
}
close $fhIn;
close $tmpFh;

my $bedCom = "windowBed ";
if ($distance)
{
        $bedCom .= "-l $distance -r $distance "
}

$bedCom = $bedCom." -a ".$tmpFile." -b ".$gff;

my $bedComOutput = `$bedCom`;
chomp $bedComOutput;


my @lines = split /\n/, $bedComOutput;

while (@lines)
{
	my $current = shift @lines;
	my @data = split /\t/, $current;
	#print "@data","\n";
	#exit;
	
	my $startDist = $data[1] - $data[8];
	my $stopDist = $data[1] - $data[9];
	
	my $dist2Gene = 0;        
   
	if ($startDist < 0)
	{
			#SNP before gene
			$dist2Gene = $startDist;
	}
	elsif ($stopDist > 0)
	{
		   $dist2Gene = $stopDist;
	}
	
	my %infos;
	my @qualif=split /;/, $data[13];
	foreach (@qualif)
	{
		my ($key,$value) = split /=/,$_;
		$value =~ s/%20/_/g;
		$value =~ s/%2/_/g;
		$infos{$key}=$value;
	}
        
	my $outline = $data[0]."\t".$data[1]."\t".$data[3]."\t".$data[4]."\t".$infos{ID}."\t".$dist2Gene."\t".$infos{Note};
	
	if ($vcf ne "null")
	{
		$outline .="\t".$vcfInfos{$data[0]}{$data[1]}{'REF'}."\t".$vcfInfos{$data[0]}{$data[1]}{'ALT'};
		$outline .="\t".$vcfInfos{$data[0]}{$data[1]}{'QUAL'}."\t".$vcfInfos{$data[0]}{$data[1]}{'ANN'};
	}
	else
	{
		$outline .="\tNA\tNA\tNA\tNA";
	}
	if ($knet ne "null")
	{
		$outline .= "\t".$knetInfos{$infos{ID}}{"SCORE"};
		if (defined $knetInfos{$infos{ID}}{"QTL"})
		{
			my $QTL = join(",",@{$knetInfos{$infos{ID}}{"QTL"}});
			$outline .= "\t$QTL";
		}
		else
		{
			$outline .="\tNA";
		}
		
		my @listEvidence=("TO","Phenotype","BioProc","Publication");
		foreach my $evidence (@listEvidence)
		{			
			if (exists $knetInfos{$infos{ID}}{"EVIDENCE"}{$evidence})
			{
					my $TO = join(",",@{$knetInfos{$infos{ID}}{"EVIDENCE"}{$evidence}});
					$outline .= "\t$TO";
			}
			else
			{
				$outline .= "\tNA";
			}
		}

	}
	else
	{
		$outline .= "\tNA\tNA\tNA\tNA\tNA\tNA";
	}
	print $fhOut $outline;
	print $fhOut "\n";
}
close $fhOut;



#Cleaning data

system("rm -Rf $tmpFile");

exit;



############################
#
# SUBS
#
############################

sub vcfReader
{
	my $headerLine = `grep -m 1 "#CHROM" $vcf`;
	chomp $headerLine;
	my @headers = split /\t/, $headerLine;
	my $nbHeader = scalar @headers;
	
	open (my $fhVCF, "<", $vcf) or (warn ("\nCannot open the $vcf file, cannot annotate the SNP...\n$!\n") && return 0);
	
	while (my $line = <$fhVCF>)
	{
		next if $line =~ m/^#/;
		next if $line =~ m/^$/;
		chomp $line;
		my @fields = split /\t/, $line;
		#Parsing the line
		my $i = 0;
		while ($i < $nbHeader)
		{
			my $localHeader = $headers[$i];
			#print $localHeader,"\n";
			#Parsing infos fields
			if ($localHeader eq "INFOS")
			{
				my @infos = split /;/, $fields[$i];
				while (@infos)
				{
					my ($key,$value) = split /=/, shift @infos;
					$vcfInfos{$fields[0]}{$fields[1]}{$key}=$value;
				}
			}
			else
			{
				$vcfInfos{$fields[0]}{$fields[1]}{$localHeader} = $fields[$i]
			}
			$i++;
		}
	}
}

sub knetReader
{
	my $conversion = &convertName;
	
	open (my $fhKnet, "<", $knet) or warn  (warn ("\nCannot open the $knet file, cannot add knet infos...\n$!\n") && return 0);
	my $j = 0;
	my @headers;
	while (my $line = <$fhKnet>)
	{
		chomp $line;
		next if $line =~ m/^#/;
		next if $line =~ m/^$/;

		if ($j == 0)
		{
			@headers = split /\t/, $line;
			$j=1;
			next;
		}
		my %tempHash;
		my @fields = split/\t/, $line;
		my $i = 0;
		my $geneID = $fields[1];
		
		if (defined $conversion->{$geneID})
		{
			$geneID = $conversion->{$geneID};
			
		}
		while (@fields)
		{
			my $local = shift @fields;
			my $type = $headers[$i];
			#print $type,"-->",$local,"\n";
			
			if ($type eq "QTL")
			{
				my @listQTL = split /\|\|/, $local;
				$tempHash{$type} = \@listQTL;
				
			}
			elsif ($type eq "EVIDENCE")
			{
				my @evidences = split /\|\|/, $local;
				while (@evidences)
				{
					my @listEv = split /\/\//, shift @evidences;
					my $subtype = shift @listEv;
					$tempHash{$type}{$subtype}=\@listEv;
				}
			}
			else
			{
				$tempHash{$type} = $local;
			}
			$i++;
		}
		if (scalar keys %{$conversion})
		{
			if (defined $$conversion{$geneID})
			{
				$geneID = $$conversion{$geneID}
			}
			
		}
		$knetInfos{$geneID}=\%tempHash;
		#last;
	}
}

sub convertName
{ #Will re-encode the knet data from the correct gff position
	open (my $fhConvert, "<", $converter) or die ("\nCannot open the converter file $converter:\n$!\n");
	my %hash;
	while (my $line = <$fhConvert>)
	{
		chomp $line;
		next if $line =~ m/^$ /;
		my @fields = split/\s+/, $line;
		my @secondFields = split /,/, $fields[1];
		my $targetOk = uc($fields[0]);
		while (@secondFields)
		{
			my $local = shift @secondFields;
			if ($local eq "None")
			{
				$local = $targetOk;
			}
			
			$local =~ s/\.\d{1,2}$//;
			$hash{$targetOk} = $local;
			
		}
	}
	return \%hash;
}