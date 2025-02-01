#TODO: configure the vizualization to follow the center of mass

proc uniqueList {list} {
set new {}
foreach item $list {
if {[lsearch $new $item] < 0} {
lappend new $item
}
}
return $new
}


proc groupContinuousIntegers {inputList} {
set result {}
set currentGroup {}
set previousValue -1

foreach value [lsort -integer $inputList] {
if {$value == $previousValue + 1} {
lappend currentGroup $value
} else {
if {[llength $currentGroup] > 0} {
lappend result $currentGroup
}
set currentGroup $value
}
set previousValue $value
}

if {[llength $currentGroup] > 0} {
lappend result $currentGroup
}

return $result
}

proc removeOverlap {list1 list2} {
set result {}
set set2 [list {*}$list2]

foreach item $list1 {
if {$item ni $set2} {
lappend result $item
}
}

return $result
}

if { [llength $argv] != 3 || [lindex $argv 0] == "-h" } {
    puts "Usage: vmd -e script.tcl -args <proteinFile> <topology> <trajectoryFile>"
    puts "Description: This script requires 3 arguments:"
    puts "  <proteinFile>   - Path to the protein structure file."
    puts "  <topologyFile>  - Path to protein topology file "
    puts "  <trajectoryFile> - Path to the trajectory file."
    puts "  <outputFile> - Path to Output file name "
}

set proteinFile [lindex $argv 0]

set topologyFile [lindex $argv 1]
set trajectoryFile [lindex $argv 2]
set outputFile [lindex $argv 3]

# Extracting relevant components from the trajectory file name
set components [split $trajectoryFile "."]
set proteinPrefix [join [lrange $components 0 1] "."]

# Loading protein structure and PSF (Protein Structure File)
mol new $proteinFile
#mol new $workingDir/$psfFile

# Loading the trajectory data into VMD
set id [mol new $topologyFile]
mol addfile $trajectoryFile waitfor all molid $id
set n_frames [molinfo $id get numframes]
animate goto 1


# Creating snapshots of the protein structure over time
set rangeMin 1
set rangeMax $n_frames

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#All of your code for decorating the protein goes here
#Decorating the protein
if {0} {
mol modselect 0 2 protein
mol modstyle 0 2 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor 0 2 ColorID 8
mol modmaterial 0 2 AOChalky
mol delete 0
mol delete 1
display ambientocclusion on
}

mol delete 0

mol modselect 0 1 protein
display resetview
mol modmaterial 0 1 AOChalky
display ambientocclusion on
mol modstyle 0 1 NewCartoon 0.300000 10.000000 4.100000 0

animate goto 1


#I'm not able to work out why my selection is not working
set helixSel [atomselect top "helix"]
set helixResNum [uniqueList [$helixSel get "residue"]]

set discreteHelixResidues [groupContinuousIntegers $helixResNum]

set betaSel [atomselect top "sheet"]
set betaResNum [uniqueList [$betaSel get "residue"]]
set discreteBetaResidues [groupContinuousIntegers $betaResNum]

puts "////////////////////////////////////////////"
puts $discreteHelixResidues
puts $discreteBetaResidues
puts "////////////////////////////////////////////"


set proteinSel [atomselect top "protein"]
set proteinResNum [uniqueList [$proteinSel get "residue"]]

set strandProtein [removeOverlap $proteinResNum $helixResNum]
set strandProtein [removeOverlap $strandProtein $betaResNum]
set strandProtein [groupContinuousIntegers $strandProtein]


set modselectionNumber 1
set colorIDNum 0
set testStructure [lindex $discreteHelixResidues 0]
puts $testStructure



#Checking if the file needs deleted
set colorStructureFileName "$outputFile"

if {[file exists $colorStructureFileName ]} {
	if {[file delete $colorStructureFileName]} {
		puts "File $colorStructureFileName  deleted successfully."
	} else {
		puts "Failed to delete $colorStructureFileName."
	}
} else {
       	puts "File $colorStructureFileName does not exist."
}


#I need this to go before the helix loop
set file_handle [open "$colorStructureFileName"  w]

puts "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
puts $discreteHelixResidues
puts "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"

#I execute the routine for coloring the secondary structures here

puts "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

foreach helix $discreteHelixResidues {
set colorIDNum [expr $colorIDNum + 1]
if {$colorIDNum == 8 || $colorIDNum == 6 || $colorIDNum == 2} {
set colorIDNum [expr $colorIDNum + 1]
}

set colorIDandStructureIndex [list $colorIDNum $helix]
puts $file_handle $colorIDandStructureIndex


puts $colorIDandStructureIndex

foreach residueNum $helix {
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID $colorIDNum
mol addrep 1
mol modselect $modselectionNumber 1 resid $residueNum
mol selection resid $residueNum
mol material Opaque
set modselectionNumber [expr $modselectionNumber + 1]
}
}


foreach sheet $discreteBetaResidues {
set colorIDNum [expr $colorIDNum + 1]
if {$colorIDNum == 8 || $colorIDNum == 6 || $colorIDNum == 2} {
set colorIDNum [expr $colorIDNum + 1]
}



#Generate the data structure here
#I just need the color ID and the indices

puts $colorIDandStructureIndex

set colorIDandStructureIndex [list $colorIDNum $sheet]
puts $file_handle $colorIDandStructureIndex

foreach residueNum $sheet {
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID $colorIDNum
mol addrep 1
mol modselect $modselectionNumber 1 resid $residueNum
mol selection resid $residueNum
mol material Opaque
set modselectionNumber [expr $modselectionNumber + 1]
}
}
puts "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"


set colorIDNum 8

mol delrep 0 1



foreach strand $strandProtein  {
foreach residueNum $strand {
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID $colorIDNum
mol addrep 1
mol modselect $modselectionNumber 1 resid $residueNum
mol selection resid $residueNum
mol material Opaque
set modselectionNumber [expr $modselectionNumber + 1]
}
}



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



#Moving to the beginning to get helical data

#Setting the initial configuration
rotate y by 90
animate goto 1
#display resetview
axes location Off



set currentProtein [atomselect top "protein"]
set COM [measure center $currentProtein weight mass]
set highestZ [lindex $COM 2]
set translationDistance 0

for {set i $rangeMin} {$i <= 13} {incr i} {
    
    if {$i % 800 == 0} {
        scale by 0.833000
    }

    animate goto $i
        
    #Collecting the center of mass
    #The command you are probably looking for is  moveby $movesel {x y z}
    set currentProtein [atomselect top "protein"]
    set COM [measure center $currentProtein weight mass]
    set z [lindex $COM 2]

    if {$z > $highestZ} {
        #Note: translation works as expected, still working this one out
         
        set cameraMovement [expr $z - $highestZ]       
        set translationDistance [expr $translationDistance + $cameraMovement]
        set negitiveTranslationDistance [expr $translationDistance * -1 ]
        set transformationArray [list 0 0 $negitiveTranslationDistance]
        $currentProtein moveby $transformationArray
        set highestZ $z

        
        } else {
             set negitiveTranslationDistance [expr $translationDistance * -1 ]
             set transformationArray [list 0 0 $negitiveTranslationDistance]
             $currentProtein moveby $transformationArray
             
        }
        # Rendering the image using TachyonLOptiXInternal renderer
        #render TachyonLOptiXInternal "$saveDir/$proteinPrefix.frame.$i.ppm" display %s

}

#exit


#TODO: generate a data structure which will allow you to determine the indices and colors
#TODO: of your protein
	#Done
