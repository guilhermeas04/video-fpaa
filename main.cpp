#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Simplex{
    private:
    int qtdLinhas;    // Quantidade de restrições (linhas da tabela Simplex)
    int qtdColunas;    // Quantidade de variáveis (colunas da tabela Simplex)
    vector<vector<float>> A;    // Matriz que armazena os coeficientes das restrições
    vector<float> B;    // Vetor que armazena os termos independentes (lado direito das restrições)
    vector<float> C;    // Vetor que armazena os coeficientes da função objetivo
    float valorMaximo;    // Armazena o valor máximo encontrado para a função objetivo
    bool isIlimitado;    // Indica se o problema é ilimitado (sem solução ótima finita)

public:
    // Construtor da classe Simplex
    // Recebe:
    // - matrix: matriz dos coeficientes das restrições
    // - b: vetor dos termos independentes
    // - c: vetor dos coeficientes da função objetivo
    Simplex(vector<vector<float>> matrix, vector<float> b, vector<float> c){
        // Inicializa o valor máximo com zero
        valorMaximo = 0;

        // Inicialmente considera que o problema não é ilimitado
        isIlimitado = false;

        // Obtém a quantidade de linhas e colunas da matriz
        qtdLinhas = matrix.size();
        qtdColunas = matrix[0].size();

        // Redimensiona as estruturas internas
        A.resize(qtdLinhas, vector<float>(qtdColunas, 0));
        B.resize(b.size());
        C.resize(c.size());

        // Copia os valores da matriz de restrições recebida
        for (int i = 0; i < qtdLinhas; i++){
            for (int j = 0; j < qtdColunas; j++){
                A[i][j] = matrix[i][j];
            }
        }

        // Copia os coeficientes da função objetivo
        for (int i = 0; i < c.size(); i++){
            C[i] = c[i];
        }

        // Copia os termos independentes das restrições
        for (int i = 0; i < b.size(); i++){
            B[i] = b[i];
        }
    }
    
    // Executa uma iteração do algoritmo Simplex.
    // Retorna true quando a solução ótima é encontrada
    // ou quando o problema é identificado como ilimitado.
    bool calculaAlgoritmoSimplex(){
        // Verifica se a tabela atual já representa uma solução ótima.
        // Caso seja ótima, não há necessidade de continuar as iterações.
        if (verificaOtimalidade() == true)
        {
            return true;
        }

        // Encontra a coluna pivô.
        // A coluna pivô corresponde ao coeficiente mais negativo
        // da função objetivo.
        int colunaPivo = findcolunaPivo();

        // Verifica se o problema foi identificado como ilimitado.
        // Nesse caso não existe solução ótima finita.
        if (isIlimitado == true)
        {
            cout << "Erro: problema ilimitado." << endl;
            return true;
        }

        // Encontra a linha pivô utilizando o teste da razão mínima.
        int linhaPivo = findlinhaPivo(colunaPivo);

        // Realiza o pivoteamento para gerar a próxima tabela Simplex.
        doPivotting(linhaPivo, colunaPivo);

        // Indica que ainda existem iterações a serem executadas.
        return false;
    }

    // Verifica se a solução atual é ótima.
    // Uma solução é considerada ótima quando todos os coeficientes
    // da função objetivo são maiores ou iguais a zero.
    bool verificaOtimalidade(){
        // Indica se a solução atual é ótima.
        bool isOptimal = false;

        // Contador de coeficientes não negativos.
        int positveValueCount = 0;

        // Percorre todos os coeficientes da função objetivo.
        for (int i = 0; i < C.size(); i++)
        {
            float value = C[i];

            // Conta quantos coeficientes são maiores ou iguais a zero.
            if (value >= 0)
            {
                positveValueCount++;
            }
        }

        // Se todos os coeficientes forem não negativos,
        // a solução encontrada é ótima.
        if (positveValueCount == C.size())
        {
            isOptimal = true;

            // Exibe a tabela final.
            print();
        }

        return isOptimal;
    }

    // Realiza a operação de pivoteamento da tabela Simplex.
// O pivoteamento transforma a coluna pivô em uma coluna identidade,
// permitindo que uma nova variável entre na base.
    void doPivotting(int linhaPivo, int colunaPivo){
        // Obtém o valor pivô localizado na interseção
        // entre a linha pivô e a coluna pivô.
        float valorPivo = A[linhaPivo][colunaPivo];

        // Armazena os valores da linha pivô original.
        float linhaPivoVals[qtdColunas];

        // Armazena os valores da coluna pivô.
        float colunaPivoVals[qtdLinhas];

        // Armazena a nova linha pivô após a normalização.
        float novaLinha[qtdColunas];

        // Atualiza o valor da função objetivo.
        valorMaximo =
            valorMaximo -
            (C[colunaPivo] * (B[linhaPivo] / valorPivo));

        // Copia os valores da linha pivô.
        for (int coluna = 0; coluna < qtdColunas; coluna++)
        {
            linhaPivoVals[coluna] = A[linhaPivo][coluna];
        }

        // Copia os valores da coluna pivô.
        for (int linha = 0; linha < qtdLinhas; linha++)
        {
            colunaPivoVals[linha] = A[linha][colunaPivo];
        }

        // Normaliza a linha pivô dividindo todos os seus
        // elementos pelo valor pivô.
        for (int coluna = 0; coluna < qtdColunas; coluna++)
        {
            novaLinha[coluna] = linhaPivoVals[coluna] / valorPivo;
        }

        // Atualiza o termo independente da linha pivô.
        B[linhaPivo] = B[linhaPivo] / valorPivo;

        // Atualiza as demais linhas da matriz de restrições.
        for (int linha = 0; linha < qtdLinhas; linha++)
        {
            // A linha pivô já foi tratada anteriormente.
            if (linha != linhaPivo)
            {
                for (int coluna = 0; coluna < qtdColunas; coluna++)
                {
                    float fatorMultiplicacao = colunaPivoVals[linha];

                    A[linha][coluna] =
                        A[linha][coluna] -
                        (fatorMultiplicacao * novaLinha[coluna]);
                }
            }
        }

        // Atualiza os termos independentes das restrições.
        for (int linha = 0; linha < B.size(); linha++)
        {
            if (linha != linhaPivo)
            {
                float fatorMultiplicacao = colunaPivoVals[linha];

                B[linha] =
                    B[linha] -
                    (fatorMultiplicacao * B[linhaPivo]);
            }
        }

        // Obtém o coeficiente da função objetivo
        // correspondente à coluna pivô.
        float fatorMultiplicacao = C[colunaPivo];

        // Atualiza a linha da função objetivo.
        for (int coluna = 0; coluna < C.size(); coluna++)
        {
            C[coluna] =
                C[coluna] -
                (fatorMultiplicacao * novaLinha[coluna]);
        }

        // Substitui a linha pivô antiga pela nova linha normalizada.
        for (int coluna = 0; coluna < qtdColunas; coluna++)
        {
            A[linhaPivo][coluna] = novaLinha[coluna];
        }
    } 
    
    // Exibe a matriz atual de restrições da tabela Simplex.
    // Utilizado para acompanhar a evolução das iterações do algoritmo.
    void print(){
        // Percorre todas as linhas da matriz.
        for (int linha = 0; linha < qtdLinhas; linha++)
        {
            // Percorre todas as colunas da matriz.
            for (int coluna = 0; coluna < qtdColunas; coluna++)
            {
                cout << A[linha][coluna] << " ";
            }

            // Quebra de linha ao final de cada linha da matriz.
            cout << endl;
        }

        // Linha em branco para melhorar a visualização da saída.
        cout << endl;
    }

    // Encontra a coluna pivô da tabela Simplex.
    // A coluna pivô é definida como a coluna que possui o menor
    // coeficiente na função objetivo. Enquanto existir um valor
    // negativo na função objetivo, ainda há possibilidade de
    // melhorar a solução atual.
    int findcolunaPivo(){
        // Inicialmente considera a primeira coluna como candidata.
        int posicao = 0;

        // Armazena o menor coeficiente encontrado na função objetivo.
        float menorValor = C[0];

        // Percorre todos os coeficientes da função objetivo.
        for (int coluna = 1; coluna < C.size(); coluna++)
        {
            // Se encontrar um valor menor, atualiza a posição da coluna pivô.
            if (C[coluna] < menorValor)
            {
                menorValor = C[coluna];
                posicao = coluna;
            }
        }

        // Retorna o índice da coluna pivô.
        return posicao;
    }

    // Encontra a linha pivô utilizando o teste da razão mínima.
    // A linha pivô é determinada dividindo os termos independentes
    // pelos coeficientes positivos da coluna pivô. A menor razão
    // positiva indica qual variável deixará a base.
    int findlinhaPivo(int colunaPivo){
        // Armazena os coeficientes positivos da coluna pivô.
        float valoresPositivos[qtdLinhas];

        // Armazena os resultados do teste da razão mínima.
        vector<float> resultado(qtdLinhas, 0);

        // Conta quantos coeficientes da coluna pivô são menores ou iguais a zero.
        int quantidadeValoresNaoPositivos = 0;

        // Identifica os coeficientes positivos da coluna pivô.
        for (int linha = 0; linha < qtdLinhas; linha++)
        {
            if (A[linha][colunaPivo] > 0)
            {
                valoresPositivos[linha] = A[linha][colunaPivo];
            }
            else
            {
                valoresPositivos[linha] = 0;
                quantidadeValoresNaoPositivos++;
            }
        }

        // Se todos os coeficientes forem não positivos,
        // o problema é ilimitado.
        if (quantidadeValoresNaoPositivos == qtdLinhas)
        {
            isIlimitado = true;
        }
        else
        {
            // Calcula as razões entre os termos independentes
            // e os coeficientes positivos da coluna pivô.
            for (int linha = 0; linha < qtdLinhas; linha++)
            {
                float coeficiente = valoresPositivos[linha];

                if (coeficiente > 0)
                {
                    resultado[linha] = B[linha] / coeficiente;
                }
                else
                {
                    resultado[linha] = 0;
                }
            }
        }

        // Procura a menor razão positiva.
        float menorRazao = 99999999;
        int indiceLinhaPivo = 0;

        for (int linha = 0; linha < resultado.size(); linha++)
        {
            if (resultado[linha] > 0)
            {
                if (resultado[linha] < menorRazao)
                {
                    menorRazao = resultado[linha];
                    indiceLinhaPivo = linha;
                }
            }
        }

        // Retorna o índice da linha pivô.
        return indiceLinhaPivo;
    }
    
    // Executa o algoritmo Simplex até encontrar a solução ótima.
    void CalculateSimplex(){
        // Controla o término das iterações.
        bool fim = false;

        // Exibe a tabela inicial.
        cout << "Tabela inicial (Nao otima)" << endl;
        print();

        cout << endl;
        cout << "Tabela final (Solucao otima)" << endl;

        // Executa as iterações do algoritmo enquanto
        // a solução ótima não for encontrada.
        while (!fim)
        {
            bool resultado = calculaAlgoritmoSimplex();

            if (resultado == true)
            {
                fim = true;
            }
        }

        // Exibe os valores das variáveis da solução ótima.
        cout << "Valores das variaveis de decisao" << endl;

        for (int coluna = 0; coluna < A.size(); coluna++)
        {
            // Conta quantos zeros existem na coluna atual.
            int quantidadeZeros = 0;

            // Armazena a linha onde está o valor 1 da coluna básica.
            int indiceLinhaBase = 0;

            for (int linha = 0; linha < qtdLinhas; linha++)
            {
                if (A[linha][coluna] == 0.0)
                {
                    quantidadeZeros++;
                }
                else if (A[linha][coluna] == 1)
                {
                    indiceLinhaBase = linha;
                }
            }

            // Se a coluna for básica (um 1 e os demais valores 0),
            // o valor da variável é obtido no vetor B.
            if (quantidadeZeros == qtdLinhas - 1)
            {
                cout << "Variavel "
                    << coluna + 1
                    << ": "
                    << B[indiceLinhaBase]
                    << endl;
            }
            else
            {
                // Variáveis não básicas possuem valor zero.
                cout << "Variavel "
                    << coluna + 1
                    << ": 0"
                    << endl;
            }
        }

        cout << endl;

        // Exibe o valor ótimo da função objetivo.
        cout << "Valor maximo da funcao objetivo: "
            << valorMaximo
            << endl;
    }
};

int main(){
    // Quantidade de colunas da matriz de restrições
    // (3 variáveis de decisão + 3 variáveis de folga)
    int quantidadeColunas = 6;

    // Quantidade de linhas da matriz de restrições
    // (número de restrições do problema)
    int quantidadeLinhas = 3;

    // Coeficientes da função objetivo
    // Max Z = 7x1 + 5x2 + 4x3
    float funcaoObjetivo[] =
    {
        -7, -5, -4, 0, 0, 0
    };

    // Termos independentes das restrições
    float termosIndependentes[] =
    {
        100,
        180,
        200
    };

    // Matriz dos coeficientes das restrições
    // Restrições:
    // x1 + x2 + x3 <= 100
    // 2x1 + x2 + 3x3 <= 180
    // x1 + 4x2 + 2x3 <= 200
    float matrizOriginal[3][6] =
    {
        {1, 1, 1, 1, 0, 0},
        {2, 1, 3, 0, 1, 0},
        {1, 4, 2, 0, 0, 1}
    };

    // Cria a matriz dinâmica que será utilizada pela classe Simplex
    vector<vector<float>> matrizRestricoes(
        quantidadeLinhas,
        vector<float>(quantidadeColunas, 0)
    );

    // Vetor dos termos independentes
    vector<float> vetorTermosIndependentes(
        quantidadeLinhas,
        0
    );

    // Vetor dos coeficientes da função objetivo
    vector<float> vetorFuncaoObjetivo(
        quantidadeColunas,
        0
    );

    // Copia os valores da matriz estática para a matriz dinâmica
    for (int linha = 0; linha < quantidadeLinhas; linha++)
    {
        for (int coluna = 0; coluna < quantidadeColunas; coluna++)
        {
            matrizRestricoes[linha][coluna] =
                matrizOriginal[linha][coluna];
        }
    }

    // Copia os termos independentes
    for (int linha = 0; linha < quantidadeLinhas; linha++)
    {
        vetorTermosIndependentes[linha] =
            termosIndependentes[linha];
    }

    // Copia os coeficientes da função objetivo
    for (int coluna = 0; coluna < quantidadeColunas; coluna++)
    {
        vetorFuncaoObjetivo[coluna] =
            funcaoObjetivo[coluna];
    }

    // Cria o objeto responsável pela execução do método Simplex
    Simplex simplex(
        matrizRestricoes,
        vetorTermosIndependentes,
        vetorFuncaoObjetivo
    );

    simplex.CalculateSimplex();

    return 0;
}
