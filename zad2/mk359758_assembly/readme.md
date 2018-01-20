# Algorytm asemblacji odczytów z sekwencjonowania

### Uruchamianie
1. `python assembly.py` `<input.fasta>` `<output.fasta>`

   * `<input.fasta>` - plik zawierający odczyty z sekwencjonowania
   * `<output.fasta>` - plik zawierający utworzone contigi

### Algorytm
1. Wykonujemy korekcję błędów dla odczytów za pomocą algorytmu opisanego na wykładzie.
    Wykorzystaliśmy implementację podaną na ćwiczeniach, z progiem równym 1:
    http://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_ErrorCorrect.ipynb

2. Budujemy graf De Bruijna dla k = 17. Do budowy grafu wykorzystaliśmy część implementacji podanej na ćwiczeniach:
    https://gist.github.com/BenLangmead/5298132

3. Zbudowany graf optymalizujemy zachłannie. Dla każdego wierzchołka
wybieramy sąsiada, który ma największą wagę (najwięcej krawędzi w multigrafie).
Każdy wierzchołek, w zoptymalizowanym grafie, będzie miał stopień wyjściowy
równy maksymalnie 1.

4. Przechodzimy graf BFS-em dopóki wszystkie wierzchołki nie będą odwiedzone. Pojedyncze wykonanie BFS-a znajduje nam ścieżkę, która staje się contigiem.

5. Eksportujemy znalezione contigi do pliku w formacie fasta.

### Parametry
1. Przetestowaliśmy różne parametry progu błędu oraz k. Oto najlepsze rezultaty
jakie udało się nam uzyskać:

| Plik     |   K  | Próg | Rezultat |
|:--------:|:----:|:----:|:--------:|
| reads1 | 17 | 2  | 0.7547 |
| reads2 | 17 | 1  | 0.2719 |
| reads3 | 14 | 1  | 0.1415 |
