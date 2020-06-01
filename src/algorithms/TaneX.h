// Strutovsky, 29.08.2019

#pragma once

#include <string>

#include "model/RelationData.h"
#include "parser/CSVParser.h"
#include "util/PositionListIndex.h"


class Tane{
private:
  CSVParser inputGenerator;
public:
  constexpr static char INPUT_FILE_CONFIG_KEY[] = "inputFile";

 //TODO: these consts should go in class (or struct) Configuration
  const double maxFdError = 0.001;
  const double maxUccError = 0.001;
  const int maxArity = 20;
  //TODO: DO NOT FORGET ABOUT MAXARITY!! - SET ARITY IN CONSTRUCTOR
  // no evidence that there will be multiple Tane instances + tricky static stuff => just const

  int countOfFD = 0;
  int countOfUCC = 0;
  long aprioriMillis = 0;
  explicit Tane(fs::path& path) : inputGenerator(path) {}
  explicit Tane(fs::path&& path) : inputGenerator(path) {} //std::move(path)
  long execute();

  static double calculateZeroAryFdError(shared_ptr<ColumnData> rhs, shared_ptr<ColumnLayoutRelationData> relationData);
  static double calculateFdError(shared_ptr<PositionListIndex> lhsPli, shared_ptr<PositionListIndex> jointPli, shared_ptr<ColumnLayoutRelationData> relationData);
  static double calculateUccError(shared_ptr<PositionListIndex> pli, shared_ptr<ColumnLayoutRelationData> relationData);
  //static double round(double error) { return ((int)(error * 32768) + 1)/ 32768.0; }
  //can't say if this method is necessary -- TODO: use this function?
  void registerFD(Vertical& lhs, shared_ptr<Column> rhs, double error, shared_ptr<RelationalSchema> schema);
  void registerFD(shared_ptr<Vertical> lhs, shared_ptr<Column> rhs, double error, shared_ptr<RelationalSchema> schema);
  void registerUCC(Vertical& key, double error, shared_ptr<RelationalSchema> schema);

};
