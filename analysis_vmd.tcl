proc avg_pos {molid sel outputfile} {
    set num_steps [molinfo $molid get numframes]
    set selection [atomselect top "[$sel text]"]
    set pos [measure avpos $sel]
    $sel set {x y z} $pos ;
    $sel writepdb $outputfile.pdb
    }

proc align_fixedatoms {molid1 molid2 sel1 sel2} {
    set num_steps [molinfo $molid2 get numframes]
    # use frame 0 for the reference
    # the frame being compared
    set compare [atomselect $molid2 "[$sel2 text]"]
    set all [atomselect $molid2 all]
    set reference [atomselect $molid1 "[$sel1 text]" frame first]

    for {set frame 0} {$frame <= $num_steps} {incr frame} {
	# get the correct frame

	$compare frame $frame
	$all frame $frame

	set trans_mat [measure fit $compare $reference]
	# do the alignment for all atoms
	$all move $trans_mat
    }
    puts "#Aligned [$sel2 text] based on fixed atoms [$sel1 text] "
}

# ============================== Measuring RMSD=========================
proc rmsd {sel molid1 molid2 name} {
    set n [molinfo $molid2 get numframes ]
    set outfile [open "$name" w ]
    #printing input given by user
    puts $outfile "#Given INPUT: \n #==================\n"
    puts $outfile "#Sel : [ $sel text ] \n #Molid : $molid1\n\n #Molid : $molid1\n\n"
    puts $outfile "#Total Frames found : $n \n\n"
    puts $outfile "#   Frame    RMSD\n==================\n"
    puts "#Given INPUT: \n==================\n"
    puts "#Sel : [$sel  text]\n # Molid : $molid1\n\n #Molid : $molid1\n\n"
    puts "#Total Frames found : $n \n\n"
    set sum 0
    #=========CYCLE STARTS==============================

    set selA [atomselect $molid1 "[$sel text]" frame 0]
    set selB [atomselect $molid2 "[$sel text]"]
    for {set i 0} {$i <= $n} {incr i} {
	$selB frame $i
	set transform_matrx [ measure fit $selB $selA ]
	#set mov_sel [atomselect 0 "all" frame $i ]
	#$mov_sel move $transform_matrx 
	$selB  move $transform_matrx
	#weighted rmsd
	set rmsdAB [measure rmsd $selA $selB weight mass ]
	#set rmsdAB [measure rmsd $selA $selB ]
	#puts [format "Frame: %5d     RMSD : %7.2f" $i $rmsdAB]
	puts $outfile [format " %5d    %7.2f" $i $rmsdAB]
	set sum [expr $sum+$rmsdAB]
    }
    set avg_rmsd [expr $sum/$n]

    # Standard Deviation
    set sum 0
    set selA [atomselect $molid1 "[$sel text]" frame 0]
    set selB [atomselect $molid2 "[$sel text]"]
    for {set i 0} {$i <= $n} {incr i} {
	$selB frame $i
	#weighted rmsd
	set rmsdAB [measure rmsd $selA $selB weight mass ]
	set Numerator [expr pow (($rmsdAB-$avg_rmsd),2)]
	set sum [expr $sum + $Numerator ]
	# sum = sum of (x-x_avg)**2
    }
    set std [expr sqrt ($sum/($n-1))]
    #puts " Standard Deviation : $std "
    #------------------------------------------------#
    puts [format "\n\n #Avg. Rmsd is        : %7.2f\n " $avg_rmsd ]
    puts [format " #Standard deviation  : %7.2f\n\n " $std ]
    puts " Data has stored in $name file . "
    puts $outfile [format "\n\n#Avg. Rmsd is : %7.2f " $avg_rmsd ]
    puts $outfile [format "\n\n#Std. is : %7.2f " $std ]
    close $outfile
    puts "Done."
}

#=============Calculating Distance====================================
proc distance {sel1 sel2 molid distance_plot histogram_plot} {
    set outfile [open "$distance_plot" w ]
    set hist_file [open "$histogram_plot" w ]
    set nbin 101
    set n [molinfo $molid get numframes]
    #--------writing given inputs--------------
    puts $outfile "\n #GIVEN INPUT DATA :\n#====================\n\n"
    puts $outfile "# Sel1 : [ $sel1 text] \n# Sel2 : [ $sel2 text ] \n\n"
    puts $outfile "# No. of Frames found : $n \n\n"
    puts  "\nGIVEN INPUT DATA :\n====================\n"
    puts  " Sel1 : [$sel1 text] \n Sel2 : [$sel2 text]\n\n"
    puts  " No. of Frames found : $n \n\n"
    set sum 0
    set atom1 [atomselect $molid "[$sel1 text]" ]
    set atom2 [atomselect $molid "[$sel2 text]" ]
    for  {set i 0}  {$i <= $n}  {incr i}  {

	$atom1 frame $i 
	set atom1_posi [ lindex [ $atom1 get "x y z"] 0 ]
	#       puts "Atom1 Positions : $atom1_posi\n"
	#	$atom1 delete

	$atom2 frame $i
	set atom2_posi [ lindex [ $atom2 get "x y z"]  0 ]
	#       puts "Atom2 Positions : $atom2_posi \n"
	#	$atom2 delete

	set dist [ vecdist $atom1_posi $atom2_posi ]
	set distance($i.r) $dist
	puts $outfile [format " %5d    %7.2f" $i  $dist]
	#       puts  "Frame : $i\t Distance : $dist "
	set sum [expr $sum + $dist]
    }
    #--------Averaging
    set Avg_dist [expr $sum/$n]
    #--------Standard Deviation
    set sum 0
    for  {set i 0}  {$i <= $n}  {incr i}  {
	set atom1 [atomselect $molid "[$sel1 text]" frame $i ]
	set atom1_posi [ lindex [ $atom1 get "x y z"] 0 ]
	$atom1 delete
	set atom2 [atomselect $molid "[$sel2 text]" frame $i ]
	set atom2_posi [ lindex [ $atom2 get "x y z"]  0 ]
	$atom2 delete
	set dist [ vecdist $atom1_posi $atom2_posi ]
	set Numarator [expr pow (($dist-$Avg_dist),2)]
	set sum [expr $sum + $Numarator ]
	# sum = sum of (x-x_avg)**2
    }
    set std [expr sqrt ($sum/($n-1))]
    #puts " Standard Deviation : $std "
    #------------------------------------------------#
    #HISTOGRAM
    set d_min distance(0.r)
    set d_max distance(0.r)
    for {set i 0} {$i < $n} {incr i} {
	if {$distance($i.r) < $d_min} {set d_min $distance($i.r)}
	if {$distance($i.r) > $d_min} {set d_max $distance($i.r)}
    }
    set width [expr ($d_max - $d_min) / ($nbin -1)]
    for {set k 0} {$k < $nbin} {incr k} {
	set distribution($k) 0
    }
    for {set i 0} {$i < $n} {incr i} {
	set k [expr int(($distance($i.r) - $d_min) / $width)]
	incr distribution($k)
    }
    for {set k 0} {$k < $nbin} {incr k} {
	puts $hist_file "[expr $d_min + $k * $width] $distribution($k)"
    }
    close  $hist_file
    #--------------------------------------------------#
    puts $outfile [ format "\n\n #Avg. Distance       : %7.2f\n\n"  $Avg_dist ]
    puts $outfile [ format "#Standard Deviation  : %7.2f\n\n"  $std ]
    puts  [ format "\n Avg. Distance       : %7.2f (A)\n"  $Avg_dist ]
    puts  [ format " Standard Deviation  : %7.2f\n\n"  $std ]
    puts  " Note : Data has been stored in files $distance_plot & $histogram_plot \n\n"
    close $outfile
}

# calculat gofr 
proc cal_gofr {sel1 sel2 delta rmax first last outputfile} {  

    set outfile [open $outputfile w]
    set gr [measure gofr $sel1 $sel2 delta $delta rmax $rmax usepbc 1 selupdate 1 first $first last $last step 1]

    set r [lindex $gr 0]
    set gr2 [lindex $gr 1]
    set igr [lindex $gr 2]
    
    set i 0
    foreach j $r k $gr2 l $igr {
	puts $outfile "$j $k $l"
    }
    puts $outfile "# Sel1 : [ $sel1 text] \n #Sel2 : [ $sel2 text ] \n\n"    
    close $outfile

}

#=============Calculating Angle====================================
proc cal_angle {sel1 sel2 sel3 molid angle_plot} {
    set outfile [open "$angle_plot" w ]
    set nbin 101
    set n [molinfo $molid get numframes]
    #--------writing given inputs--------------
    puts $outfile "\n #GIVEN INPUT DATA :\n#====================\n\n"
    puts $outfile "# Sel1 : [ $sel1 text] \n# Sel2 : [ $sel2 text ] \n Sel3 : [$sel3 text]\n\n"
    puts $outfile "# No. of Frames found : $n \n\n"
    puts  "\nGIVEN INPUT DATA :\n====================\n"
    puts  " Sel1 : [$sel1 text] \n Sel2 : [$sel2 text] \n \n\n"
    puts  " No. of Frames found : $n \n\n"
    set sum 0
    set atom1 [atomselect $molid "[$sel1 text]" ]
    set atom2 [atomselect $molid "[$sel2 text]" ]
    set atom3 [atomselect $molid "[$sel3 text]" ]
    for  {set i 0}  {$i <= $n}  {incr i}  {
	$atom1 frame $i
	$atom2 frame $i
	$atom3 frame $i
        set A1 [ $sel1 list ]
    	set A2 [ $sel2 list ]
        set A3 [ $sel3 list ]
	set list [ subst {$A1 $A2 $A3} ]
	set ang [ measure angle $list frame $i]
	set cal_angle($i.r) $ang
	puts $outfile [format " %5d    %7.2f" $i  $ang]
    }
    close $outfile
}

#=============Calculating dihedralangle====================================
proc cal_dihedral {sel1 sel2 sel3 sel4 molid angle_plot} {
    set outfile [open "$angle_plot" w ]
    set nbin 101
    set n [molinfo $molid get numframes]
    #--------writing given inputs--------------
    puts $outfile "\n #GIVEN INPUT DATA :\n#====================\n\n"
    puts $outfile "# Sel1 : [ $sel1 text] \n# Sel2 : [ $sel2 text ] \n# Sel3 : [$sel3 text] \n #Sel4 : [$sel4 text] \n\n"
    puts $outfile "# No. of Frames found : $n \n\n"
    puts  "\nGIVEN INPUT DATA :\n====================\n"
    puts  " Sel1 : [$sel1 text] \n Sel2 : [$sel2 text] \n Sel3 : [$sel3 text] \n Sel4 : [$sel4 text]\n\n"
    puts  " No. of Frames found : $n \n\n"
    set sum 0
    set atom1 [atomselect $molid "[$sel1 text]" ]
    set atom2 [atomselect $molid "[$sel2 text]" ]
    set atom3 [atomselect $molid "[$sel3 text]" ]
    set atom4 [atomselect $molid "[$sel4 text]" ]
    for  {set i 0}  {$i <= $n}  {incr i}  {
	$atom1 frame $i
	$atom2 frame $i
	$atom3 frame $i
	$atom4 frame $i
	set A1 [ $sel1 list ]
    	set A2 [ $sel2 list ]
        set A3 [ $sel3 list ]
        set A4 [ $sel4 list ]
	set list [ subst {$A1 $A2 $A3 $A4} ]
	set dihed [ measure dihed $list frame $i]
	set cal_dihedral($i.r) $dihed
	puts $outfile [format " %5d  %7.2f" $i  $dihed]
    }
    close $outfile
}

proc cal_sasa {solute output} {
    set num_steps [molinfo top get numframes]
    set fout [open $output w]
    set nf [molinfo top get numframes] 
    set solute [atomselect top "[ $solute text ]"]
    set all [atomselect top all]		

    for {set frame 0} {$frame <= $nf} {incr frame} {

    	$solute frame $frame
	$all frame $frame
	set sasavalue [measure sasa 1.40 $solute] 
	puts $fout "$frame $sasavalue"
    }	
}

proc cal_sasarestrict {solute ressel output} {
    set num_steps [molinfo top get numframes]
    set fout [open $output w]
    set nf [molinfo top get numframes] 
    set solute [atomselect top "[ $solute text ]"]
    set ressel [atomselect top "[ $ressel text ]"]	
    set all [atomselect top all]		

    for {set frame 0} {$frame <= $nf} {incr frame} {

    	$solute frame $frame
	$ressel frame $frame
	$all frame $frame
	set sasavalue [measure sasa 1.40 $solute -restrict $ressel] 
	puts $fout "$frame $sasavalue"
    }	
}
