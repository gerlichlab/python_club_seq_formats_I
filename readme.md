# Python club about sequencing formats I

In this first python club about sequencing formats we will talk about formats that store nucleotide level information. These file are often text-based (due to easy of early implementation and later convention) and have a wide range of information content.

## FASTA files

The FASTA format is as easy as it gets. It consists of one or multiple read entries, each of which conforms to this format:

```
    header line ---- >READNAME ADDITIONALINFORMATION
    Sequence    ---- ACGTCCGTGAGAGAGTCGTGCACGCGACGT
```
The only restriction that it enforces is that the `READNAME` - the first set of characters after the `>` sign needs to be unique. But of course this restriction is delegated to the programs that work with FASTA files and they might or might not tell you that this happened.

FASTA files are often used in situations where quality information does not make a lot of sense, for example when describing reference genomes or primers or plasmids.

## Running code in this problem set

In order work with the material here, you need to install [docker](https://docs.docker.com/get-docker/). You can then start a notebook server with all the necessary libraries installed with `docker-compose up` from the command line in the repo directory. Then just open up a browser and navigate to http://localhost:9999. This will prompt you for a password, which will be disclosed in the lecture. You can then run code and edit files in your browser.

Alternatively, you can open any of the notebooks in vscode, open up the command palette and specify a remote notebook, by pasting `http://localhost:9999` as follows:

* Open notebook file
* Hit `ctrl+shift+P` or `magic-apple-button + shift+P` to open the command palette
    ![asdf](pictures\command_palette.png)
* Type in `Specify local or remote jupyter server` 
    ![asdf](pictures\specify_local_or_remote.png)
* Click on `Existing`
    ![asdf](pictures\existing.png)
* Type in `http://localhost:9999` and hit enter.
    ![asdf](pictures\localhost.png)

If you want to run the tests for the three functions that you need to implement, run the following command form within the repo directory `docker-compose -f run_tests.yml up`.