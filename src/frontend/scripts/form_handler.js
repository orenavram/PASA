
function openAdvanced() {
    var flag = document.getElementById("open_advanced_flag");
    if (flag.value === "off"){
        flag.value = "on"
        var x = document.getElementById("enrichment_threshold_div");
        x.style.display = "block";
        if (document.getElementById("example_page").value === "no") {
            var checkbox = document.getElementById("checkbox");
            if (checkbox.checked == false) {
                // no enzyme
                x = document.getElementById("enzyme_div");
                x.style.display = "block";
            }
        } else {
                x = document.getElementById("enzyme_div");
                x.style.display = "block";
        }
    } else {
        flag.value = "off"
        // show enzyme
        var x = document.getElementById("enzyme_div");
        x.style.display = "none";
        x = document.getElementById("enrichment_threshold_div");
        x.style.display = "none";
    }
}


function RawOrTxt() {
  var peptides_txt = document.getElementById("checkbox");
  // If the peptides_txt checkbox is checked, display the peptides text options
  if (peptides_txt.checked == true){
      // show peptides lists input
    document.getElementById("el_peptides_div").style.display = "block";
    document.getElementById("ft_peptides_div").style.display = "block";
    document.getElementById("el_raw_div").style.display = "none";
    document.getElementById("ft_raw_div").style.display = "none";
    document.getElementById("enzyme_div").style.display = "none";
    document.getElementById("el").value = '';
    document.getElementById("ft").value = '';
    //   if (document.getElementById("enrichment_threshold_div").style.display === "none"){
    // }
  } else {
      // show raw input
    document.getElementById("el_peptides_div").style.display = "none";
    document.getElementById("ft_peptides_div").style.display = "none";
    document.getElementById("el_raw_div").style.display = "block";
    document.getElementById("ft_raw_div").style.display = "block";
    if (document.getElementById("open_advanced_flag").value == "on") {
        document.getElementById("enzyme_div").style.display = "block";
        document.getElementById("enzyme").value = 'Trypsin';
    }
    document.getElementById("el_peptides").value = '';
    document.getElementById("ft_peptides").value = '';
  }
}




function validate_read_files() {
    var raws = [document.getElementsByName("elution_data_1")[0].value,
        document.getElementsByName("elution_data_2")[0].value,
        document.getElementsByName("elution_data_3")[0].value,
        document.getElementsByName("flowthrough_data_1")[0].value,
        document.getElementsByName("flowthrough_data_2")[0].value,
        document.getElementsByName("flowthrough_data_3")[0].value]

    var db = document.getElementsByName("database")[0].value

    // for (i = 0; i < raws.length; i++){
    //     if (!(raws[i].endsWith('raw'))) {
    //         alert("One of the raw files is illegal. Only raw formats are allowed.");
    //         return false;
    //     }
    // }

    if (!(db.endsWith('fasta') || db.endsWith('gz') || db.endsWith('zip') || db.endsWith('rar') || db.endsWith('.tar.gz'))) {
        alert("DB file is illegal. Only zip/tar.gz/fasta/rar/gz formats are allowed.");
        return false;
    }

    return true;

}





