
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
    if (!document.getElementsByName("db").value){
        alert("DB files ZIP is missing!");
        return false;
    }
    var peptides_txt = document.getElementById("checkbox");
    if (peptides_txt.checked == true){
        // peptides lists input
        if (!document.getElementsByName("el_peptides_div").value){
            alert("Elution peptides file is missing!");
            return false;
        }
        if (!document.getElementsByName("ft_peptides_div").value){
            alert("Flow-through peptides file is missing!");
            return false;
        }
    } else {
      // raw input
        if (!document.getElementsByName("el_raw_div").value){
            alert("Elution raw files ZIP is missing!");
            return false;
        }
        if (!document.getElementsByName("ft_raw_div").value){
            alert("Flow-through raw files ZIP is missing!");
            return false;
        }
    }
    return true;
}





