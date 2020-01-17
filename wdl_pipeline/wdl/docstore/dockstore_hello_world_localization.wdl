version 1.0


workflow HelloWorldLocalization {
    input {
        File FILE
        Array[File] FILES
    }

    call head as single_h {
        input:
            myFile=FILE
    }

    scatter (file in FILES) {
        call head as scatter_h {
            input:
                myFile=file
        }
    }

    output {
        File fileOut = scatter_h.myHead
        Array[File] filesOut = scatter_h.myHead
    }
}

task head {
    input {
        File myFile
    }
    command <<<
        head ~{myFile}
        ln -s ~{myFile} mySymlinkFile
        head mySymlinkFile >>output
    >>>
    output {
        File myHead = read_file("output")
    }
    runtime {
        docker: "tpesout/vgp_minimap2:latest"
        cpu: 1
    }

}

