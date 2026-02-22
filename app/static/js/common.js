// Shared utility functions

function setError(m) {
    document.getElementById("error-text").innerHTML = m;
}

// Collapsible toggle setup
function initCollapsibles() {
    var coll = document.getElementsByClassName("collapsible");
    for (var i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function() {
            this.classList.toggle("active-collapsible");
            var content = this.nextElementSibling;
            if (content.style.maxHeight) {
                content.style.maxHeight = null;
            } else {
                content.style.maxHeight = content.scrollHeight + "px";
            }
        });
    }
}

document.addEventListener("DOMContentLoaded", function() {
    initCollapsibles();
});
