while(my $line=<>)
{
	my @cur=split/\s+/, $line;
	print $#cur."\n";
}
