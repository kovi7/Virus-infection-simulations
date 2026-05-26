# Projekt: przestrzenne modelowanie infekcji wirusowej w Wolfram Mathematica

## Autorzy
- Agata Paluch
- Basia Pawłowska
- Justyna Kowalska

## 1. Cel projektu

Główna idea projektu:

1. Najpierw pokazujemy prosty automat komórkowy, w którym infekcja przechodzi tylko przez zakażonych sąsiadów.
2. Potem rozszerzamy model o wolne wiriony `V(t,x,y)`, dyfuzję, fazę opóźnienia produkcji wirusa, chemotaksję, przeciwciała, leczenie, szczepienie i mutację.
3. Na końcu porównujemy, jak różne założenia matematyczne zmieniają przebieg symulacji.

## 2. Zawartość paczki

W paczce znajdują się:

- `FinalVirusModel.wl`, główny plik z kodem modelu i wszystkimi funkcjami eksperymentów.
- `01_Model_podstawowy_i_typy_wirusow.nb`, pierwszy notebook, model prosty, porównanie modelu prostego z zaawansowanym, profile fenotypowe.
- `02_Dyfuzja_transmisja_lag.nb`, drugi notebook, dyfuzja wirionów, transmisja `cell-free` i `cell-to-cell`, faza `lag`.
- `03_Odpornosc_chemotaksja_pamiec.nb`, trzeci notebook, układ odpornościowy, ruch losowy vs chemotaksja, szybka i opóźniona odpowiedź, pamięć immunologiczna.
- `04_Terapia_szczepienie_mutacje_Hill.nb`, czwarty notebook, leczenie, szczepienie, strategie dawkowania, mutacja ucieczkowa, heatmapa dla 10 symulacji i funkcja Hilla.

## 3. Jak uruchomić projekt

1. Rozpakuj ZIP do jednego folderu.
2. Otwórz wybrany notebook `.nb` w Wolfram Mathematica.
3. Na początku notebooka uruchom komórkę:

```wolfram
Get[FileNameJoin[{NotebookDirectory[], "FinalVirusModel.wl"}]];
```

4. Uruchamiaj komórki po kolei od góry do dołu.
5. Jeżeli Mathematica pamięta stare definicje funkcji, wykonaj:

```wolfram
Quit[]
```

Następnie uruchom notebook od początku.

## 4. Najważniejsze funkcje

### Model prosty

```wolfram
RunSimpleScenario[name]
```

Uruchamia prosty automat komórkowy. Przykład:

```wolfram
simple = RunSimpleScenario["SimpleBaseline"];
AnimateResult[simple]
PlotCounts[simple]
```

### Model zaawansowany

```wolfram
RunScenario[name]
```

Uruchamia model hybrydowy, czyli automat komórkowy plus pola `V`, `C`, `A`. Przykład:

```wolfram
adv = RunScenario["BaselineAdvanced"];
AnimateResult[adv, "ShowVirus" -> True]
CompactSummaryPlot[adv]
```

### Gotowe eksperymenty

W projekcie przygotowane są gotowe funkcje:

```wolfram
ExperimentSimpleVsAdvanced[]
ExperimentVirusTypes[]
ExperimentDiffusion[]
ExperimentTransmissionModes[]
ExperimentLagPhase[]
ExperimentChemotaxis[]
ExperimentImmuneProfiles[]
ExperimentMemory[]
ExperimentDrugVaccination[]
ExperimentDosing[]
ExperimentMutation[]
ExperimentTherapyHeatmap[]
ExperimentHill[]
```

Każda funkcja zwraca `Association`, czyli strukturę z wynikami, animacjami i wykresami.

Przykład:

```wolfram
exp = ExperimentTransmissionModes[];
exp["Animation"]
exp["InfectedPlot"]
exp["VirusPlot"]
```

## 5. Ważne założenie: start infekcji w centrum

W aktualnej wersji projektu animowane symulacje startują od centrum siatki.

Jeżeli `InitialInfected -> 1`, zakażona jest komórka centralna. Jeżeli `InitialInfected > 1`, jedna komórka startuje w centrum, a pozostałe są rozmieszczane losowo.

To jest warunek początkowy wybrany dla czytelności. Nie oznacza, że biologicznie wirus zawsze zaczyna w środku tkanki. Taki start pozwala łatwo zobaczyć falę infekcji, wpływ dyfuzji i różnice między scenariuszami.

## 6. Legenda kolorów

W animacjach używana jest następująca konwencja:

- biały, komórka zdrowa i podatna na zakażenie,
- pomarańczowy, komórka zakażona w fazie `lag`, jeszcze bez intensywnej produkcji wirionów,
- czerwony, komórka zakażona aktywnie, produkująca wiriony,
- czarny, komórka martwa po infekcji,
- jasnozielony, komórka chroniona, na przykład przez szczepienie,
- żółty, mutant ucieczkowy,
- fioletowy, komórka usunięta przez leczenie albo toksycznie uszkodzona przez terapię,
- fioletowe punkty, wolne wiriony próbkowane z pola `V(t,x,y)`,
- niebieskie punkty, komórki odpornościowe,
- niebieskawe tło, obecność leku w scenariuszach terapii.

Kolory są konwencją wizualną, nie rzeczywistymi barwami biologicznych struktur.

## 7. Model prosty a model zaawansowany

### Model prosty

W modelu prostym infekcja zależy tylko od liczby zakażonych sąsiadów. Dla zdrowej komórki liczymy `nInf`, czyli liczbę zakażonych sąsiadów. Prawdopodobieństwo zakażenia wynosi:

```text
P = 1 - (1 - BetaContact)^nInf
```

Model prosty jest czytelny, ale pomija wolne wiriony, dyfuzję, opóźnienie produkcji wirusa i odpowiedź immunologiczną.

### Model zaawansowany

Model zaawansowany jest hybrydowy. Łączy:

- dyskretne stany komórek,
- ciągłe pola na siatce.

Najważniejsze pola:

```text
V(t,x,y), wolne wiriony,
C(t,x,y), sygnał chemotaktyczny,
A(t,x,y), przeciwciała.
```

W modelu zaawansowanym zakażenie zachodzi dwiema drogami:

1. `cell-to-cell`, przez kontakt z zakażonym sąsiadem,
2. `cell-free`, przez wolne wiriony dyfundujące po siatce.

## 8. Leczenie w aktualnej wersji

Lek działa w trzech mechanizmach:

1. `DrugProductionEffect`, zmniejsza produkcję nowych wirionów przez zakażone komórki.
2. `DrugInfectionEffect`, zmniejsza efektywną zakaźność w transmisji `cell-free` i `cell-to-cell`.
3. `DrugInfectedKilling`, może bezpośrednio usuwać zakażone komórki, które wtedy są oznaczane jako fioletowe.

Dodatkowo `DrugToxicity` może oznaczać uproszczony efekt uboczny, czyli uszkodzenie zdrowych komórek przez terapię.

Niebieskawe tło w animacji nie jest osobnym polem PDE. To wizualizacja aktualnego poziomu leku `drug` w symulacji.

## 9. Scenariusze w notebookach

### Notebook 01

Pokazuje:

- prosty model,
- model zaawansowany,
- porównanie obu modeli przy podobnej lokalnej zakaźności,
- profile fenotypowe, czyli szybki wirus o wysokiej transmisji vs wolniejszy i bardziej uszkadzający.


### Notebook 02

Pokazuje:

- wpływ dyfuzji wirionów,
- różnicę między transmisją `cell-free` i `cell-to-cell`,
- fazę `lag`, czyli opóźnienie produkcji wirionów,
- widok komórki plus wolne wiriony, inspirowany Figure 2a z pracy Graw i Perelson.

### Notebook 03

Pokazuje:

- pole chemotaktyczne `C(t,x,y)`,
- ruch losowy i chemotaktyczny komórek odpornościowych,
- szybką/silną i opóźnioną/słabszą odpowiedź odpornościową,
- pamięć immunologiczną, modelowaną jako szybszy start odpowiedzi i początkowy poziom przeciwciał.

### Notebook 04

Pokazuje:

- leczenie,
- szczepienie,
- wczesną i spóźnioną interwencję,
- dawkę jednorazową i dawkowanie cykliczne,
- mutację ucieczkową,
- heatmapę dla 10 symulacji
- funkcję Hilla
- barierę krew-mózg

## 10. Ograniczenia projektu

Najważniejsze ograniczenia:

- model jest dwuwymiarowy,
- parametry są przeskalowane do symulacji dydaktycznej,
- nie modelujemy całego organizmu,
- nie modelujemy układu krwionośnego, temperatury, wielu narządów ani pełnej immunologii, jedynie barierę krew-mózg,
- `SARS-CoV-2-like` i `Hantavirus-like` to profile modelowe, nie rzeczywiste kliniczne symulacje chorób,
- wyniki należy interpretować jakościowo, a nie ilościowo.

## 11. Literatura i inspiracje

Projekt opiera się koncepcyjnie na dwóch pracach:

1. Graw F., Perelson A. S., `Modeling Viral Spread`, Annual Review of Virology, 2016.
2. Cai Y., Zhao Z., Zhuge C., `The spatial dynamics of immune response upon virus infection through hybrid dynamical computational model`, Frontiers in Immunology, 2023.

Pierwsza praca uzasadnia rozdzielenie transmisji `cell-free` i `cell-to-cell`. Druga praca uzasadnia podejście hybrydowe, czyli automat komórkowy połączony z polami reakcji-dyfuzji.
